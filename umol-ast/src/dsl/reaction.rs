//! Reaction DSL.
//!
//! `ReactionDsl` wraps a `ReactionAst` together with the `ReactionMetadata` that records
//! the surface-form bindings (the lhs molecule metadata plus created-entity keywords and
//! atom-alias bindings). The EDN form is a map keyed by `:lhs` (a molecule map) and
//! `:deltas` (a vector of `:add` / `:remove` / `:modify` / `:constraint` operations). Each
//! entity delegates to its own entity DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use bimap::BiBTreeMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_perm::Permutation;

use super::aromatic::{AromaticSystemDsl, AromaticSystemUpdateDsl};
use super::atom::{lower_atom, raise_atom, AtomDsl, AtomUpdateDsl};
use super::bond::{lower_bond, raise_bond, BondDsl, BondUpdateDsl};
use super::config::{DeltaDefaults, ReactionDefaults};
use super::constraint::{read_constraint_dsl, ConstraintDsl};
use super::dative::{DativeBondDsl, DativeBondUpdateDsl};
use super::edn_utils::{
    consume_single_key_map_close, missing, parse_single_key_map, parse_vec,
    read_single_key_map_header, read_vec, single_key_map,
};
use super::error::ParseError;
use super::metadata::{MoleculeMetadata, ReactionMetadata};
use super::molecule::{
    parse_aromatic_system_entry, parse_atom_aliases, parse_atom_entry, parse_bond_entry,
    parse_dative_bond_entry, parse_molecule_input, parse_multicenter_bond_entry,
    parse_noncovalent_bond_entry, parse_stereo_atom_entry, parse_stereo_bond_entry,
    read_aromatic_system_entry, read_atom_aliases, read_atom_entry, read_bond_entry,
    read_dative_bond_entry, read_molecule_input, read_multicenter_bond_entry,
    read_noncovalent_bond_entry, read_stereo_atom_entry, read_stereo_bond_entry,
    render_aromatic_entry, render_dative_entry, render_molecule_edn, render_multicenter_entry,
    render_noncovalent_entry, render_stereo_atom_entry, render_stereo_bond_entry,
    render_stereo_ligand, resolve_atom_spec, AromaticSystemEntryInput, AtomEntryInput,
    BondEntryInput, DativeBondEntryInput, MoleculeDsl, MoleculeInput, MulticenterBondEntryInput,
    NoncovalentBondEntryInput, StereoAtomEntryInput, StereoBondEntryInput,
};
use super::multicenter::{MulticenterBondDsl, MulticenterBondUpdateDsl};
use super::namespace::{MoleculeNamespace, Namespace};
use super::noncovalent::{NoncovalentBondDsl, NoncovalentBondUpdateDsl};
use super::refs::{
    read_aromatic_system_ref, read_atom_ref, read_bond_ref, read_dative_bond_ref,
    read_multicenter_bond_ref, read_noncovalent_bond_ref, read_stereo_atom_ref,
    read_stereo_bond_ref, AromaticSystemRef, AtomRef, BondRef, DativeBondRef, MulticenterBondRef,
    NoncovalentBondRef, StereoAtomRef, StereoBondRef,
};
use super::stereo::{
    parse_permutation, render_edn_stereo_kind, stereo_kind_from_name, StereoAtomDsl,
    StereoAtomUpdateDsl, StereoBondDsl, StereoBondUpdateDsl,
};
use crate::ast::atom::{AtomAst, AtomUpdate};
use crate::ast::bond::{BondAst, BondUpdate};
use crate::ast::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, ConstraintDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta,
};
use crate::ast::edit::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use crate::ast::entity::Entity;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::reaction::ReactionAst;
use crate::ast::stereo::{StereoConfigurationAst, StereoConfigurationUpdate};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::{
    AromaticSystemUpdate, DativeBondUpdate, MulticenterBondUpdate, NoncovalentBondUpdate,
    StereoAtomUpdate, StereoBondUpdate, StereoKind, StereoLigand,
};

/// Surface DSL for a reaction. Pairs `ReactionAst` with `ReactionMetadata`; fields are
/// private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionDsl {
    ast: ReactionAst,
    metadata: ReactionMetadata,
}

impl ReactionDsl {
    pub fn from_parts(ast: ReactionAst, metadata: ReactionMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &ReactionAst {
        &self.ast
    }

    pub fn metadata(&self) -> &ReactionMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (ReactionAst, ReactionMetadata) {
        (self.ast, self.metadata)
    }
}

impl FromAst<ReactionAst> for ReactionDsl {
    type Ctx = ReactionDefaults;

    fn from_ast(ast: &ReactionAst, cfg: &Self::Ctx) -> Self {
        let lhs = MoleculeDsl::from_ast(&ast.lhs, &cfg.molecule_defaults())
            .into_parts()
            .0;
        let delta_cfg = cfg.delta_defaults();
        let mut deltas = ast.deltas.clone();
        for delta in deltas.iter_mut() {
            lower_delta(delta, &delta_cfg);
        }
        ReactionDsl {
            ast: ReactionAst { lhs, deltas },
            metadata: ReactionMetadata::default(),
        }
    }
}

impl IntoAst<ReactionAst> for ReactionDsl {
    type Ctx = ReactionDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> ReactionAst {
        let ReactionAst { lhs, mut deltas } = self.ast;
        let lhs = MoleculeDsl::from_parts(lhs, MoleculeMetadata::default())
            .into_ast(&cfg.molecule_defaults());
        let delta_cfg = cfg.delta_defaults();
        for delta in deltas.iter_mut() {
            raise_delta(delta, &delta_cfg);
        }
        ReactionAst { lhs, deltas }
    }
}

/// Lower a delta's embedded entity AST to DSL-display form (AST → DSL).
fn lower_delta(delta: &mut Delta, cfg: &DeltaDefaults) {
    match delta {
        Delta::Atom(AtomDelta::Add { ast, .. } | AtomDelta::Remove { ast, .. }) => {
            lower_atom(ast, &cfg.atom)
        }
        Delta::Bond(BondDelta::Add { ast, .. } | BondDelta::Remove { ast, .. }) => {
            lower_bond(ast, &cfg.bond)
        }
        Delta::DativeBond(
            DativeBondDelta::Add { ast, .. } | DativeBondDelta::Remove { ast, .. },
        ) => *ast = DativeBondDsl::from_ast(ast, &cfg.dative_bond).0,
        Delta::AromaticSystem(
            AromaticSystemDelta::Add { ast, .. } | AromaticSystemDelta::Remove { ast, .. },
        ) => *ast = AromaticSystemDsl::from_ast(ast, &cfg.aromatic_system).0,
        Delta::MulticenterBond(
            MulticenterBondDelta::Add { ast, .. } | MulticenterBondDelta::Remove { ast, .. },
        ) => *ast = MulticenterBondDsl::from_ast(ast, &cfg.multicenter_bond).0,
        Delta::NoncovalentBond(
            NoncovalentBondDelta::Add { ast, .. } | NoncovalentBondDelta::Remove { ast, .. },
        ) => *ast = NoncovalentBondDsl::from_ast(ast, &cfg.noncovalent_bond).0,
        Delta::StereoAtom(
            StereoAtomDelta::Add { ast, .. } | StereoAtomDelta::Remove { ast, .. },
        ) => *ast = StereoAtomDsl::from_ast(ast, &cfg.stereo_atom).0,
        Delta::StereoBond(
            StereoBondDelta::Add { ast, .. } | StereoBondDelta::Remove { ast, .. },
        ) => *ast = StereoBondDsl::from_ast(ast, &cfg.stereo_bond).0,
        _ => {}
    }
}

/// Raise a delta's embedded entity AST from DSL-display form (DSL → AST).
fn raise_delta(delta: &mut Delta, cfg: &DeltaDefaults) {
    match delta {
        Delta::Atom(AtomDelta::Add { ast, .. } | AtomDelta::Remove { ast, .. }) => {
            raise_atom(ast, &cfg.atom)
        }
        Delta::Bond(BondDelta::Add { ast, .. } | BondDelta::Remove { ast, .. }) => {
            raise_bond(ast, &cfg.bond)
        }
        Delta::DativeBond(
            DativeBondDelta::Add { ast, .. } | DativeBondDelta::Remove { ast, .. },
        ) => *ast = DativeBondDsl(ast.clone()).into_ast(&cfg.dative_bond),
        Delta::AromaticSystem(
            AromaticSystemDelta::Add { ast, .. } | AromaticSystemDelta::Remove { ast, .. },
        ) => *ast = AromaticSystemDsl(ast.clone()).into_ast(&cfg.aromatic_system),
        Delta::MulticenterBond(
            MulticenterBondDelta::Add { ast, .. } | MulticenterBondDelta::Remove { ast, .. },
        ) => *ast = MulticenterBondDsl(ast.clone()).into_ast(&cfg.multicenter_bond),
        Delta::NoncovalentBond(
            NoncovalentBondDelta::Add { ast, .. } | NoncovalentBondDelta::Remove { ast, .. },
        ) => *ast = NoncovalentBondDsl(ast.clone()).into_ast(&cfg.noncovalent_bond),
        Delta::StereoAtom(
            StereoAtomDelta::Add { ast, .. } | StereoAtomDelta::Remove { ast, .. },
        ) => *ast = StereoAtomDsl(ast.clone()).into_ast(&cfg.stereo_atom),
        Delta::StereoBond(
            StereoBondDelta::Add { ast, .. } | StereoBondDelta::Remove { ast, .. },
        ) => *ast = StereoBondDsl(ast.clone()).into_ast(&cfg.stereo_bond),
        _ => {}
    }
}

/// The reaction's resolution namespace: the lhs molecule's namespace, the delta namespace continuing
/// its id space (holding every entity a delta binds — its own alias map stays empty), and the
/// reaction's top-level atom aliases in a field of their own. A ref resolves against the union: a
/// keyword or participant key is looked up in `deltas` first, then `lhs` (the id spaces are disjoint, so at
/// most one hits). Counts come from `deltas`, which continues `lhs` and so carries the reaction-wide
/// total — the single running counter that hands out delta ids on `register_*`.
pub struct ReactionNamespace {
    lhs: MoleculeNamespace,
    deltas: MoleculeNamespace,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}

/// Where a reaction id came from: an entity of the lhs molecule, or one introduced by a delta.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EntityScope {
    Lhs,
    Deltas,
}

impl ReactionNamespace {
    fn new(lhs: MoleculeNamespace) -> Self {
        let deltas = MoleculeNamespace::continuation(&lhs);
        Self {
            lhs,
            deltas,
            atom_aliases: BiBTreeMap::new(),
        }
    }

    /// The namespace of an already-resolved reaction: the lhs molecule's namespace plus every entity
    /// an `Add` delta introduces, registered anonymously with its participants in delta order (which
    /// reproduces the per-kind delta ids). Refs resolve against it as they did at parse time.
    pub fn from_ast(reaction: &ReactionAst) -> Self {
        let free = "anonymous delta entity registration never collides";
        let mut ns = Self::new(MoleculeNamespace::from_ast(&reaction.lhs));
        for delta in reaction.deltas.iter() {
            match delta {
                Delta::Atom(AtomDelta::Add { .. }) => {
                    ns.register_atom(None).expect(free);
                }
                Delta::Bond(BondDelta::Add { atoms, .. }) => {
                    ns.register_bond(None, atoms[0], atoms[1]).expect(free);
                }
                Delta::DativeBond(DativeBondDelta::Add {
                    donors, acceptor, ..
                }) => {
                    ns.register_dative_bond(None, donors, *acceptor)
                        .expect(free);
                }
                Delta::AromaticSystem(AromaticSystemDelta::Add { atoms, .. }) => {
                    ns.register_aromatic_system(None, atoms).expect(free);
                }
                Delta::MulticenterBond(MulticenterBondDelta::Add { atoms, .. }) => {
                    ns.register_multicenter_bond(None, atoms).expect(free);
                }
                Delta::NoncovalentBond(NoncovalentBondDelta::Add { atoms, .. }) => {
                    ns.register_noncovalent_bond(None, atoms[0], atoms[1])
                        .expect(free);
                }
                Delta::StereoAtom(StereoAtomDelta::Add { site, ligands, .. }) => {
                    ns.register_stereo_atom(None, *site, ligands).expect(free);
                }
                Delta::StereoBond(StereoBondDelta::Add { site, ligands, .. }) => {
                    ns.register_stereo_bond(None, *site, ligands).expect(free);
                }
                _ => {}
            }
        }
        ns
    }

    /// Whether a keyword is free in the lhs + reaction-alias scope. The delta scope is checked by the
    /// delegated `deltas.register_*`; the two together cover the whole reaction namespace.
    fn check_keyword_free(&self, keyword: Option<&str>) -> Result<(), ParseError> {
        match keyword {
            Some(kw) if self.lhs.contains_keyword(kw) || self.atom_aliases.contains_left(kw) => {
                Err(ParseError::DuplicateKeyword(kw.to_string()))
            }
            _ => Ok(()),
        }
    }

    fn register_atom(&mut self, keyword: Option<String>) -> Result<AtomId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_atom(keyword)
    }
    fn register_bond(
        &mut self,
        keyword: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> Result<BondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_bond(keyword, a, b)
    }
    fn register_dative_bond(
        &mut self,
        keyword: Option<String>,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Result<DativeBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_dative_bond(keyword, donors, acceptor)
    }
    fn register_aromatic_system(
        &mut self,
        keyword: Option<String>,
        atoms: &[AtomId],
    ) -> Result<AromaticSystemId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_aromatic_system(keyword, atoms)
    }
    fn register_multicenter_bond(
        &mut self,
        keyword: Option<String>,
        atoms: &[AtomId],
    ) -> Result<MulticenterBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_multicenter_bond(keyword, atoms)
    }
    fn register_noncovalent_bond(
        &mut self,
        keyword: Option<String>,
        a: AtomId,
        b: AtomId,
    ) -> Result<NoncovalentBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_noncovalent_bond(keyword, a, b)
    }
    fn register_stereo_atom(
        &mut self,
        keyword: Option<String>,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Result<StereoAtomId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_stereo_atom(keyword, site, ligands)
    }
    fn register_stereo_bond(
        &mut self,
        keyword: Option<String>,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Result<StereoBondId, ParseError> {
        self.check_keyword_free(keyword.as_deref())?;
        self.deltas.register_stereo_bond(keyword, site, ligands)
    }

    /// Bind a top-level reaction atom alias, erroring if the name is already taken (any entity or
    /// alias, lhs or reaction) or the atom-spec is already aliased (bijectivity).
    fn register_atom_alias(&mut self, name: String, dsl: Box<AtomDsl>) -> Result<(), ParseError> {
        if self.contains_keyword(&name) {
            return Err(ParseError::DuplicateKeyword(name));
        }
        if self.alias_targets(&dsl) {
            return Err(ParseError::InvalidValue(
                "atom-aliases must be bijective: two names map to the same atom".into(),
            ));
        }
        self.atom_aliases.insert(name, dsl);
        Ok(())
    }

    /// Whether `dsl` is already some reaction or lhs alias's target — the bijection check.
    fn alias_targets(&self, dsl: &AtomDsl) -> bool {
        self.atom_aliases
            .iter()
            .any(|(_, existing)| existing.as_ref() == dsl)
            || self.lhs.atom_aliases().any(|(_, existing)| existing == dsl)
    }

    pub(super) fn lhs(&self) -> &MoleculeNamespace {
        &self.lhs
    }
    pub(super) fn deltas(&self) -> &MoleculeNamespace {
        &self.deltas
    }
    pub(super) fn atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(name, dsl)| (name.as_str(), dsl.as_ref()))
    }

    /// Classify a reaction id by kind: `Lhs` if its index is below the lhs count for that kind
    /// (an lhs entity), `Deltas` if at or above (a delta introduced it).
    fn scope(index: usize, lhs_count: usize) -> EntityScope {
        if index < lhs_count {
            EntityScope::Lhs
        } else {
            EntityScope::Deltas
        }
    }
    fn atom_scope(&self, id: AtomId) -> EntityScope {
        Self::scope(id.index(), self.lhs.atom_count())
    }
    fn bond_scope(&self, id: BondId) -> EntityScope {
        Self::scope(id.index(), self.lhs.bond_count())
    }
    fn dative_bond_scope(&self, id: DativeBondId) -> EntityScope {
        Self::scope(id.index(), self.lhs.dative_bond_count())
    }
    fn aromatic_system_scope(&self, id: AromaticSystemId) -> EntityScope {
        Self::scope(id.index(), self.lhs.aromatic_system_count())
    }
    fn multicenter_bond_scope(&self, id: MulticenterBondId) -> EntityScope {
        Self::scope(id.index(), self.lhs.multicenter_bond_count())
    }
    fn noncovalent_bond_scope(&self, id: NoncovalentBondId) -> EntityScope {
        Self::scope(id.index(), self.lhs.noncovalent_bond_count())
    }
    fn stereo_atom_scope(&self, id: StereoAtomId) -> EntityScope {
        Self::scope(id.index(), self.lhs.stereo_atom_count())
    }
    fn stereo_bond_scope(&self, id: StereoBondId) -> EntityScope {
        Self::scope(id.index(), self.lhs.stereo_bond_count())
    }
}

impl Namespace for ReactionNamespace {
    fn atom_count(&self) -> usize {
        self.deltas.atom_count()
    }
    fn bond_count(&self) -> usize {
        self.deltas.bond_count()
    }
    fn dative_bond_count(&self) -> usize {
        self.deltas.dative_bond_count()
    }
    fn aromatic_system_count(&self) -> usize {
        self.deltas.aromatic_system_count()
    }
    fn multicenter_bond_count(&self) -> usize {
        self.deltas.multicenter_bond_count()
    }
    fn noncovalent_bond_count(&self) -> usize {
        self.deltas.noncovalent_bond_count()
    }
    fn stereo_atom_count(&self) -> usize {
        self.deltas.stereo_atom_count()
    }
    fn stereo_bond_count(&self) -> usize {
        self.deltas.stereo_bond_count()
    }

    fn find_atom_by_keyword(&self, keyword: &str) -> Option<AtomId> {
        self.deltas
            .find_atom_by_keyword(keyword)
            .or_else(|| self.lhs.find_atom_by_keyword(keyword))
    }
    fn find_bond_by_keyword(&self, keyword: &str) -> Option<BondId> {
        self.deltas
            .find_bond_by_keyword(keyword)
            .or_else(|| self.lhs.find_bond_by_keyword(keyword))
    }
    fn find_dative_bond_by_keyword(&self, keyword: &str) -> Option<DativeBondId> {
        self.deltas
            .find_dative_bond_by_keyword(keyword)
            .or_else(|| self.lhs.find_dative_bond_by_keyword(keyword))
    }
    fn find_aromatic_system_by_keyword(&self, keyword: &str) -> Option<AromaticSystemId> {
        self.deltas
            .find_aromatic_system_by_keyword(keyword)
            .or_else(|| self.lhs.find_aromatic_system_by_keyword(keyword))
    }
    fn find_multicenter_bond_by_keyword(&self, keyword: &str) -> Option<MulticenterBondId> {
        self.deltas
            .find_multicenter_bond_by_keyword(keyword)
            .or_else(|| self.lhs.find_multicenter_bond_by_keyword(keyword))
    }
    fn find_noncovalent_bond_by_keyword(&self, keyword: &str) -> Option<NoncovalentBondId> {
        self.deltas
            .find_noncovalent_bond_by_keyword(keyword)
            .or_else(|| self.lhs.find_noncovalent_bond_by_keyword(keyword))
    }
    fn find_stereo_atom_by_keyword(&self, keyword: &str) -> Option<StereoAtomId> {
        self.deltas
            .find_stereo_atom_by_keyword(keyword)
            .or_else(|| self.lhs.find_stereo_atom_by_keyword(keyword))
    }
    fn find_stereo_bond_by_keyword(&self, keyword: &str) -> Option<StereoBondId> {
        self.deltas
            .find_stereo_bond_by_keyword(keyword)
            .or_else(|| self.lhs.find_stereo_bond_by_keyword(keyword))
    }

    fn find_bond_by_participants(&self, first: AtomId, second: AtomId) -> Option<BondId> {
        self.deltas
            .find_bond_by_participants(first, second)
            .or_else(|| self.lhs.find_bond_by_participants(first, second))
    }
    fn find_dative_bond_by_participants(
        &self,
        donors: &[AtomId],
        acceptor: AtomId,
    ) -> Option<DativeBondId> {
        self.deltas
            .find_dative_bond_by_participants(donors, acceptor)
            .or_else(|| self.lhs.find_dative_bond_by_participants(donors, acceptor))
    }
    fn find_aromatic_system_by_participants(&self, atoms: &[AtomId]) -> Option<AromaticSystemId> {
        self.deltas
            .find_aromatic_system_by_participants(atoms)
            .or_else(|| self.lhs.find_aromatic_system_by_participants(atoms))
    }
    fn find_multicenter_bond_by_participants(&self, atoms: &[AtomId]) -> Option<MulticenterBondId> {
        self.deltas
            .find_multicenter_bond_by_participants(atoms)
            .or_else(|| self.lhs.find_multicenter_bond_by_participants(atoms))
    }
    fn find_noncovalent_bond_by_participants(
        &self,
        first: AtomId,
        second: AtomId,
    ) -> Option<NoncovalentBondId> {
        self.deltas
            .find_noncovalent_bond_by_participants(first, second)
            .or_else(|| {
                self.lhs
                    .find_noncovalent_bond_by_participants(first, second)
            })
    }
    fn find_stereo_atom_by_participants(
        &self,
        site: AtomId,
        ligands: &[StereoLigand],
    ) -> Option<StereoAtomId> {
        self.deltas
            .find_stereo_atom_by_participants(site, ligands)
            .or_else(|| self.lhs.find_stereo_atom_by_participants(site, ligands))
    }
    fn find_stereo_bond_by_participants(
        &self,
        site: BondId,
        ligands: &[StereoLigand],
    ) -> Option<StereoBondId> {
        self.deltas
            .find_stereo_bond_by_participants(site, ligands)
            .or_else(|| self.lhs.find_stereo_bond_by_participants(site, ligands))
    }

    fn contains_keyword(&self, keyword: &str) -> bool {
        self.deltas.contains_keyword(keyword)
            || self.lhs.contains_keyword(keyword)
            || self.atom_aliases.contains_left(keyword)
    }

    fn find_atom_alias(&self, name: &str) -> Option<&AtomDsl> {
        self.atom_aliases
            .get_by_left(name)
            .map(|dsl| dsl.as_ref())
            .or_else(|| self.lhs.find_atom_alias(name))
    }
}

/// One unresolved delta parsed from a `:deltas` entry. Refs stay symbolic
/// (`AtomRef`/`BondRef`) and an `:add` carries the molecule entry verbatim
/// (bare atom / alias / created keyword resolved in R7); a `:modify` RHS is a
/// entity update value.
#[derive(Debug, PartialEq)]
pub(crate) enum DeltaInput {
    AtomAdd(AtomEntryInput),
    AtomRemove(AtomRef),
    AtomModify(AtomRef, AtomUpdate),
    BondAdd(BondEntryInput),
    BondRemove(BondRef),
    BondModify(BondRef, BondUpdate),
    DativeBondAdd(DativeBondEntryInput),
    DativeBondRemove(DativeBondRef),
    DativeBondModify(DativeBondRef, DativeBondUpdate),
    AromaticSystemAdd(AromaticSystemEntryInput),
    AromaticSystemRemove(AromaticSystemRef),
    AromaticSystemModify(AromaticSystemRef, AromaticSystemUpdate),
    MulticenterBondAdd(MulticenterBondEntryInput),
    MulticenterBondRemove(MulticenterBondRef),
    MulticenterBondModify(MulticenterBondRef, MulticenterBondUpdate),
    NoncovalentBondAdd(NoncovalentBondEntryInput),
    NoncovalentBondRemove(NoncovalentBondRef),
    NoncovalentBondModify(NoncovalentBondRef, NoncovalentBondUpdate),
    StereoAtomAdd(StereoAtomEntryInput),
    StereoAtomRemove(StereoAtomRef),
    StereoAtomModify(StereoAtomRef, StereoAtomUpdate),
    StereoAtomSwap(StereoAtomRef, StereoKind),
    StereoAtomMirror(StereoAtomRef, StereoKind),
    StereoAtomApply(StereoAtomRef, StereoKind, Permutation),
    StereoBondAdd(StereoBondEntryInput),
    StereoBondRemove(StereoBondRef),
    StereoBondModify(StereoBondRef, StereoBondUpdate),
    StereoBondSwap(StereoBondRef, StereoKind),
    StereoBondMirror(StereoBondRef, StereoKind),
    StereoBondApply(StereoBondRef, StereoKind, Permutation),
    ConstraintAdd(ConstraintDsl),
    ConstraintRemove(ConstraintDsl),
}

/// Raw parse target for a reaction: the lhs molecule input plus the unresolved
/// deltas. Resolution (`into_ast`, R7) lifts this to `(ReactionAst, ReactionMetadata)`.
#[derive(Debug, PartialEq)]
pub(crate) struct ReactionInput {
    lhs: MoleculeInput,
    atom_aliases: Vec<(String, Box<AtomDsl>)>,
    deltas: Vec<DeltaInput>,
}

impl ReactionInput {
    pub(crate) fn into_ast(self) -> Result<(ReactionAst, ReactionMetadata), ParseError> {
        let ReactionInput {
            lhs,
            atom_aliases,
            deltas,
        } = self;
        let (lhs, lhs_namespace) = lhs.into_ast()?;

        // The single resolution namespace: lhs entities, the delta namespace continuing its id space,
        // and the reaction's top-level aliases. Every ref — entity and constraint — resolves against
        // it; `register_*` advances the delta counter and `ReactionMetadata` is projected from it at
        // the end. No forward refs: only entities bound earlier in delta order are visible.
        let mut ns = ReactionNamespace::new(lhs_namespace);

        // Top-level reaction aliases, bijective (a name colliding with an lhs alias, or a target
        // already aliased, errors); they resolve in delta atom-specs alongside the lhs aliases.
        for (name, dsl) in atom_aliases {
            ns.register_atom_alias(name, dsl)?;
        }

        let mut resolved = Deltas::new();
        for delta in deltas {
            match delta {
                DeltaInput::AtomAdd(entry) => {
                    let ast = resolve_atom_spec(entry.spec, &ns)?;
                    let id = ns.register_atom(entry.keyword)?;
                    resolved.push(Delta::Atom(AtomDelta::Add { id, ast }));
                }
                DeltaInput::AtomRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "atom",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::Atom(AtomDelta::Remove {
                        id,
                        ast: lhs.atom(id).ast.clone(),
                    }));
                }
                DeltaInput::AtomModify(r, update) => {
                    let id = r.resolve(&ns)?;
                    if ns.atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "atom",
                            index: id.index(),
                        });
                    }
                    for d in AtomDelta::for_update(id, lhs.atom(id).ast, &update) {
                        resolved.push(Delta::Atom(d));
                    }
                }
                DeltaInput::BondAdd(entry) => {
                    let a = entry.first.resolve(&ns)?;
                    let b = entry.second.resolve(&ns)?;
                    let id = ns.register_bond(entry.keyword, a, b)?;
                    resolved.push(Delta::Bond(BondDelta::Add {
                        id,
                        atoms: [a, b],
                        ast: entry.bond.0,
                    }));
                }
                DeltaInput::BondRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "bond",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::Bond(BondDelta::Remove {
                        id,
                        atoms: lhs.bond(id).atom_ids(),
                        ast: lhs.bond(id).ast.clone(),
                    }));
                }
                DeltaInput::BondModify(r, update) => {
                    let id = r.resolve(&ns)?;
                    if ns.bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "bond",
                            index: id.index(),
                        });
                    }
                    for d in BondDelta::for_update(id, lhs.bond(id).ast, &update) {
                        resolved.push(Delta::Bond(d));
                    }
                }
                DeltaInput::DativeBondAdd(entry) => {
                    let donors = entry
                        .donors
                        .into_iter()
                        .map(|d| d.resolve(&ns))
                        .collect::<Result<Vec<_>, _>>()?;
                    let acceptor = entry.acceptor.resolve(&ns)?;
                    let id = ns.register_dative_bond(entry.keyword, &donors, acceptor)?;
                    resolved.push(Delta::DativeBond(DativeBondDelta::Add {
                        id,
                        donors,
                        acceptor,
                        ast: entry.bond.0,
                    }));
                }
                DeltaInput::DativeBondRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.dative_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "dative bond",
                            index: id.index(),
                        });
                    }
                    let view = lhs.dative_bond(id);
                    resolved.push(Delta::DativeBond(DativeBondDelta::Remove {
                        id,
                        donors: view.donor_ids().collect(),
                        acceptor: view.acceptor_id(),
                        ast: lhs.dative_bond(id).ast.clone(),
                    }));
                }
                DeltaInput::DativeBondModify(r, update) => {
                    let id = r.resolve(&ns)?;
                    if ns.dative_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "dative bond",
                            index: id.index(),
                        });
                    }
                    for d in DativeBondDelta::for_update(id, lhs.dative_bond(id).ast, &update) {
                        resolved.push(Delta::DativeBond(d));
                    }
                }
                DeltaInput::AromaticSystemAdd(entry) => {
                    let atoms = entry
                        .atoms
                        .into_iter()
                        .map(|a| a.resolve(&ns))
                        .collect::<Result<Vec<_>, _>>()?;
                    let id = ns.register_aromatic_system(entry.keyword, &atoms)?;
                    resolved.push(Delta::AromaticSystem(AromaticSystemDelta::Add {
                        id,
                        atoms,
                        ast: entry.system.0,
                    }));
                }
                DeltaInput::AromaticSystemRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.aromatic_system_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "aromatic system",
                            index: id.index(),
                        });
                    }
                    let view = lhs.aromatic_system(id);
                    resolved.push(Delta::AromaticSystem(AromaticSystemDelta::Remove {
                        id,
                        atoms: view.atom_ids().collect(),
                        ast: lhs.aromatic_system(id).ast.clone(),
                    }));
                }
                DeltaInput::AromaticSystemModify(r, rhs) => {
                    let id = r.resolve(&ns)?;
                    if ns.aromatic_system_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "aromatic system",
                            index: id.index(),
                        });
                    }
                    for d in AromaticSystemDelta::for_update(id, lhs.aromatic_system(id).ast, &rhs)
                    {
                        resolved.push(Delta::AromaticSystem(d));
                    }
                }
                DeltaInput::MulticenterBondAdd(entry) => {
                    let atoms = entry
                        .atoms
                        .into_iter()
                        .map(|a| a.resolve(&ns))
                        .collect::<Result<Vec<_>, _>>()?;
                    let id = ns.register_multicenter_bond(entry.keyword, &atoms)?;
                    resolved.push(Delta::MulticenterBond(MulticenterBondDelta::Add {
                        id,
                        atoms,
                        ast: entry.bond.0,
                    }));
                }
                DeltaInput::MulticenterBondRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.multicenter_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "multicenter bond",
                            index: id.index(),
                        });
                    }
                    let view = lhs.multicenter_bond(id);
                    resolved.push(Delta::MulticenterBond(MulticenterBondDelta::Remove {
                        id,
                        atoms: view.atom_ids().collect(),
                        ast: lhs.multicenter_bond(id).ast.clone(),
                    }));
                }
                DeltaInput::MulticenterBondModify(r, rhs) => {
                    let id = r.resolve(&ns)?;
                    if ns.multicenter_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "multicenter bond",
                            index: id.index(),
                        });
                    }
                    for d in
                        MulticenterBondDelta::for_update(id, lhs.multicenter_bond(id).ast, &rhs)
                    {
                        resolved.push(Delta::MulticenterBond(d));
                    }
                }
                DeltaInput::NoncovalentBondAdd(entry) => {
                    let first = entry.first.resolve(&ns)?;
                    let second = entry.second.resolve(&ns)?;
                    let id = ns.register_noncovalent_bond(entry.keyword, first, second)?;
                    resolved.push(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
                        id,
                        atoms: [first, second],
                        ast: entry.bond.0,
                    }));
                }
                DeltaInput::NoncovalentBondRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.noncovalent_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "noncovalent bond",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::NoncovalentBond(NoncovalentBondDelta::Remove {
                        id,
                        atoms: lhs.noncovalent_bond(id).atom_ids(),
                        ast: lhs.noncovalent_bond(id).ast.clone(),
                    }));
                }
                DeltaInput::NoncovalentBondModify(r, rhs) => {
                    let id = r.resolve(&ns)?;
                    if ns.noncovalent_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "noncovalent bond",
                            index: id.index(),
                        });
                    }
                    for d in
                        NoncovalentBondDelta::for_update(id, lhs.noncovalent_bond(id).ast, &rhs)
                    {
                        resolved.push(Delta::NoncovalentBond(d));
                    }
                }
                DeltaInput::StereoAtomAdd(entry) => {
                    let site = entry.site.resolve(&ns)?;
                    let ligands = entry
                        .ligands
                        .into_iter()
                        .map(|l| Ok(StereoLigand::new(l.atom.resolve(&ns)?, l.kind)))
                        .collect::<Result<Vec<_>, ParseError>>()?;
                    let id = ns.register_stereo_atom(entry.keyword, site, &ligands)?;
                    resolved.push(Delta::StereoAtom(StereoAtomDelta::Add {
                        id,
                        site,
                        ligands,
                        ast: entry.stereo.0,
                    }));
                }
                DeltaInput::StereoAtomRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "stereo atom",
                            index: id.index(),
                        });
                    }
                    let view = lhs.stereo_atom(id);
                    resolved.push(Delta::StereoAtom(StereoAtomDelta::Remove {
                        id,
                        site: view.site_id(),
                        ligands: view
                            .ligands()
                            .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                            .collect(),
                        ast: lhs.stereo_atom(id).ast.clone(),
                    }));
                }
                DeltaInput::StereoAtomModify(r, rhs) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "stereo atom",
                            index: id.index(),
                        });
                    }
                    for d in StereoAtomDelta::for_update(id, lhs.stereo_atom(id).ast, &rhs) {
                        resolved.push(Delta::StereoAtom(d));
                    }
                }
                DeltaInput::StereoAtomSwap(r, kind) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "transform",
                            kind: "stereo atom",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::StereoAtom(StereoAtomDelta::Swap { id, kind }));
                }
                DeltaInput::StereoAtomMirror(r, kind) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "transform",
                            kind: "stereo atom",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::StereoAtom(StereoAtomDelta::Mirror { id, kind }));
                }
                DeltaInput::StereoAtomApply(r, kind, permutation) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_atom_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "transform",
                            kind: "stereo atom",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::StereoAtom(StereoAtomDelta::Apply {
                        id,
                        kind,
                        permutation,
                    }));
                }
                DeltaInput::StereoBondAdd(entry) => {
                    let site = entry.site.resolve(&ns)?;
                    let ligands = entry
                        .ligands
                        .into_iter()
                        .map(|l| Ok(StereoLigand::new(l.atom.resolve(&ns)?, l.kind)))
                        .collect::<Result<Vec<_>, ParseError>>()?;
                    let id = ns.register_stereo_bond(entry.keyword, site, &ligands)?;
                    resolved.push(Delta::StereoBond(StereoBondDelta::Add {
                        id,
                        site,
                        ligands,
                        ast: entry.stereo.0,
                    }));
                }
                DeltaInput::StereoBondRemove(r) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "remove",
                            kind: "stereo bond",
                            index: id.index(),
                        });
                    }
                    let view = lhs.stereo_bond(id);
                    resolved.push(Delta::StereoBond(StereoBondDelta::Remove {
                        id,
                        site: view.site_id(),
                        ligands: view
                            .ligands()
                            .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                            .collect(),
                        ast: lhs.stereo_bond(id).ast.clone(),
                    }));
                }
                DeltaInput::StereoBondModify(r, rhs) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "modify",
                            kind: "stereo bond",
                            index: id.index(),
                        });
                    }
                    for d in StereoBondDelta::for_update(id, lhs.stereo_bond(id).ast, &rhs) {
                        resolved.push(Delta::StereoBond(d));
                    }
                }
                DeltaInput::StereoBondSwap(r, kind) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "transform",
                            kind: "stereo bond",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::StereoBond(StereoBondDelta::Swap { id, kind }));
                }
                DeltaInput::StereoBondMirror(r, kind) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "transform",
                            kind: "stereo bond",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::StereoBond(StereoBondDelta::Mirror { id, kind }));
                }
                DeltaInput::StereoBondApply(r, kind, permutation) => {
                    let id = r.resolve(&ns)?;
                    if ns.stereo_bond_scope(id) == EntityScope::Deltas {
                        return Err(ParseError::DeltaTargetAdded {
                            action: "transform",
                            kind: "stereo bond",
                            index: id.index(),
                        });
                    }
                    resolved.push(Delta::StereoBond(StereoBondDelta::Apply {
                        id,
                        kind,
                        permutation,
                    }));
                }
                DeltaInput::ConstraintAdd(dsl) => {
                    let c = dsl.into_ast(&ns)?;
                    resolved.push(Delta::Constraint(ConstraintDelta::Add(c)));
                }
                DeltaInput::ConstraintRemove(dsl) => {
                    let c = dsl.into_ast(&ns)?;
                    resolved.push(Delta::Constraint(ConstraintDelta::Remove(c)));
                }
            }
        }
        let metadata = ReactionMetadata::from(&ns);
        Ok((
            ReactionAst {
                lhs,
                deltas: resolved,
            },
            metadata,
        ))
    }
}

impl<'de> FromEdn<'de> for ReactionDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let input = parse_reaction_input(edn)?;
        let (ast, metadata) = input
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionDsl::from_parts(ast, metadata))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let ri = read_reaction_input(&mut de)?;
        de.expect_eof()?;
        let (ast, metadata) = ri.into_ast().map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(ReactionDsl::from_parts(ast, metadata))
    }
}

impl FromStr for ReactionDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ReactionDsl::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl ToEdn for ReactionDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_reaction_edn(&self.ast, &self.metadata)
    }
}

impl Display for ReactionDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for ReactionAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        ReactionDsl::from_edn(edn).map(|dsl| dsl.into_parts().0)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        ReactionDsl::from_edn_str(input).map(|dsl| dsl.into_parts().0)
    }
}

impl FromStr for ReactionAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for ReactionAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

/// Direct EDN rendering for `ReactionAst`. The lhs and refs emit canonical positional
/// form (no entity keywords, no aliases) since the AST carries no metadata. For keyword/alias-bearing
/// surface output, wrap in [`ReactionDsl`] with the appropriate [`ReactionMetadata`].
impl ToEdn for ReactionAst {
    fn to_edn(&self) -> Edn<'static> {
        render_reaction_edn(self, &ReactionMetadata::default())
    }
}

fn read_reaction_input(de: &mut EdnStreamDeserializer<'_>) -> Result<ReactionInput, EdnError> {
    de.consume_byte(b'{')?;
    let mut lhs = None;
    let mut atom_aliases = Vec::new();
    let mut deltas = Vec::new();
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?.into_owned();
        match key.as_str() {
            "lhs" => lhs = Some(read_molecule_input(de)?),
            "atom-aliases" => atom_aliases = read_atom_aliases(de)?,
            "deltas" => deltas = read_vec(de, read_delta_input)?,
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["reaction".into()],
                }
                .into())
            }
        }
    }
    Ok(ReactionInput {
        lhs: lhs.ok_or_else(|| missing("lhs", "reaction"))?,
        atom_aliases,
        deltas,
    })
}

fn parse_reaction_input(edn: &Edn<'_>) -> Result<ReactionInput, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "reaction map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let mut lhs = None;
    let mut atom_aliases = Vec::new();
    let mut deltas = Vec::new();
    for (k, v) in m.iter() {
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["reaction".into()],
            });
        };
        match key.name() {
            "lhs" => lhs = Some(parse_molecule_input(v)?),
            "atom-aliases" => atom_aliases = parse_atom_aliases(v)?,
            "deltas" => deltas = parse_vec(v, ":deltas", parse_delta_input)?,
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["reaction".into()],
                })
            }
        }
    }
    Ok(ReactionInput {
        lhs: lhs.ok_or(DeError::MissingField {
            key: "lhs".to_string(),
            path: vec!["reaction".into()],
        })?,
        atom_aliases,
        deltas,
    })
}

fn read_delta_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let entity = read_single_key_map_header(de)?;
    let input = match entity.as_str() {
        "atom" => read_delta_atom_input(de)?,
        "bond" => read_delta_bond_input(de)?,
        "dative-bond" => read_delta_dative_bond_input(de)?,
        "aromatic-system" => read_delta_aromatic_system_input(de)?,
        "multicenter-bond" => read_delta_multicenter_bond_input(de)?,
        "noncovalent-bond" => read_delta_noncovalent_bond_input(de)?,
        "stereo-atom" => read_delta_stereo_atom_input(de)?,
        "stereo-bond" => read_delta_stereo_bond_input(de)?,
        "constraint" => read_delta_constraint_input(de)?,
        e => return Err(DeError::Custom(format!("unknown reaction delta :{e}")).into()),
    };
    consume_single_key_map_close(de, "delta")?;
    Ok(input)
}

fn read_delta_atom_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::AtomAdd(read_atom_entry(de)?),
        "remove" => DeltaInput::AtomRemove(read_atom_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_atom_ref(de)?;
            let s = de.read_string()?;
            let dsl: AtomUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("atom-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("atom :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::AtomModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown atom delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "atom delta")?;
    Ok(input)
}

fn read_delta_bond_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::BondAdd(read_bond_entry(de)?),
        "remove" => DeltaInput::BondRemove(read_bond_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_bond_ref(de)?;
            let s = de.read_string()?;
            let dsl: BondUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("bond-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("bond :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::BondModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown bond delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "bond delta")?;
    Ok(input)
}

fn read_stereo_kind(de: &mut EdnStreamDeserializer<'_>) -> Result<StereoKind, EdnError> {
    let name = de.read_keyword_name()?;
    stereo_kind_from_name(name.as_ref()).map_err(Into::into)
}

fn read_delta_dative_bond_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::DativeBondAdd(read_dative_bond_entry(de)?),
        "remove" => DeltaInput::DativeBondRemove(read_dative_bond_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_dative_bond_ref(de)?;
            let s = de.read_string()?;
            let dsl: DativeBondUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("dative-bond-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("dative-bond :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::DativeBondModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown dative-bond delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "dative-bond delta")?;
    Ok(input)
}

fn read_delta_aromatic_system_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::AromaticSystemAdd(read_aromatic_system_entry(de)?),
        "remove" => DeltaInput::AromaticSystemRemove(read_aromatic_system_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_aromatic_system_ref(de)?;
            let s = de.read_string()?;
            let dsl: AromaticSystemUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("aromatic-system-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("aromatic-system :modify expects [ref dsl]".into()).into(),
                );
            }
            DeltaInput::AromaticSystemModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown aromatic-system delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "aromatic-system delta")?;
    Ok(input)
}

fn read_delta_multicenter_bond_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::MulticenterBondAdd(read_multicenter_bond_entry(de)?),
        "remove" => DeltaInput::MulticenterBondRemove(read_multicenter_bond_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_multicenter_bond_ref(de)?;
            let s = de.read_string()?;
            let dsl: MulticenterBondUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("multicenter-bond-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("multicenter-bond :modify expects [ref dsl]".into()).into(),
                );
            }
            DeltaInput::MulticenterBondModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown multicenter-bond delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "multicenter-bond delta")?;
    Ok(input)
}

fn read_delta_noncovalent_bond_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::NoncovalentBondAdd(read_noncovalent_bond_entry(de)?),
        "remove" => DeltaInput::NoncovalentBondRemove(read_noncovalent_bond_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_noncovalent_bond_ref(de)?;
            let s = de.read_string()?;
            let dsl: NoncovalentBondUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("noncovalent-bond-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("noncovalent-bond :modify expects [ref dsl]".into()).into(),
                );
            }
            DeltaInput::NoncovalentBondModify(r, dsl.0)
        }
        o => return Err(DeError::Custom(format!("unknown noncovalent-bond delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "noncovalent-bond delta")?;
    Ok(input)
}

fn read_delta_stereo_atom_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::StereoAtomAdd(read_stereo_atom_entry(de)?),
        "remove" => DeltaInput::StereoAtomRemove(read_stereo_atom_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_atom_ref(de)?;
            let s = de.read_string()?;
            let dsl: StereoAtomUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("stereo-atom-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("stereo-atom :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::StereoAtomModify(r, dsl.0)
        }
        "swap" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_atom_ref(de)?;
            let kind = read_stereo_kind(de)?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("stereo-atom :swap expects [ref kind]".into()).into());
            }
            DeltaInput::StereoAtomSwap(r, kind)
        }
        "mirror" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_atom_ref(de)?;
            let kind = read_stereo_kind(de)?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("stereo-atom :mirror expects [ref kind]".into()).into(),
                );
            }
            DeltaInput::StereoAtomMirror(r, kind)
        }
        "apply" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_atom_ref(de)?;
            let kind = read_stereo_kind(de)?;
            let s = de.read_string()?;
            let permutation = parse_permutation(s.as_ref(), kind.degree())
                .map_err(|e| DeError::subgrammar("permutation", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("stereo-atom :apply expects [ref kind cycles]".into()).into(),
                );
            }
            DeltaInput::StereoAtomApply(r, kind, permutation)
        }
        o => return Err(DeError::Custom(format!("unknown stereo-atom delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "stereo-atom delta")?;
    Ok(input)
}

fn read_delta_stereo_bond_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::StereoBondAdd(read_stereo_bond_entry(de)?),
        "remove" => DeltaInput::StereoBondRemove(read_stereo_bond_ref(de)?),
        "modify" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_bond_ref(de)?;
            let s = de.read_string()?;
            let dsl: StereoBondUpdateDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("stereo-bond-update", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("stereo-bond :modify expects [ref dsl]".into()).into());
            }
            DeltaInput::StereoBondModify(r, dsl.0)
        }
        "swap" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_bond_ref(de)?;
            let kind = read_stereo_kind(de)?;
            if !de.try_consume_byte(b']')? {
                return Err(DeError::Custom("stereo-bond :swap expects [ref kind]".into()).into());
            }
            DeltaInput::StereoBondSwap(r, kind)
        }
        "mirror" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_bond_ref(de)?;
            let kind = read_stereo_kind(de)?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("stereo-bond :mirror expects [ref kind]".into()).into(),
                );
            }
            DeltaInput::StereoBondMirror(r, kind)
        }
        "apply" => {
            de.consume_byte(b'[')?;
            let r = read_stereo_bond_ref(de)?;
            let kind = read_stereo_kind(de)?;
            let s = de.read_string()?;
            let permutation = parse_permutation(s.as_ref(), kind.degree())
                .map_err(|e| DeError::subgrammar("permutation", e))?;
            if !de.try_consume_byte(b']')? {
                return Err(
                    DeError::Custom("stereo-bond :apply expects [ref kind cycles]".into()).into(),
                );
            }
            DeltaInput::StereoBondApply(r, kind, permutation)
        }
        o => return Err(DeError::Custom(format!("unknown stereo-bond delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "stereo-bond delta")?;
    Ok(input)
}

fn read_delta_constraint_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let op = read_single_key_map_header(de)?;
    let input = match op.as_str() {
        "add" => DeltaInput::ConstraintAdd(read_constraint_dsl(de)?),
        "remove" => DeltaInput::ConstraintRemove(read_constraint_dsl(de)?),
        o => return Err(DeError::Custom(format!("unknown constraint delta op :{o}")).into()),
    };
    consume_single_key_map_close(de, "constraint delta")?;
    Ok(input)
}

fn parse_delta_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (entity, body) = parse_single_key_map(edn, "delta")?;
    match entity {
        "atom" => parse_delta_atom_input(body),
        "bond" => parse_delta_bond_input(body),
        "dative-bond" => parse_delta_dative_bond_input(body),
        "aromatic-system" => parse_delta_aromatic_system_input(body),
        "multicenter-bond" => parse_delta_multicenter_bond_input(body),
        "noncovalent-bond" => parse_delta_noncovalent_bond_input(body),
        "stereo-atom" => parse_delta_stereo_atom_input(body),
        "stereo-bond" => parse_delta_stereo_bond_input(body),
        "constraint" => parse_delta_constraint_input(body),
        e => Err(DeError::Custom(format!("unknown reaction delta :{e}"))),
    }
}

fn parse_delta_atom_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "atom delta")?;
    match op {
        "add" => Ok(DeltaInput::AtomAdd(parse_atom_entry(payload)?)),
        "remove" => Ok(DeltaInput::AtomRemove(AtomRef::from_edn(payload)?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "atom :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["atom delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "atom :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::AtomModify(
                AtomRef::from_edn(&v[0])?,
                AtomUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!("unknown atom delta op :{o}"))),
    }
}

fn parse_delta_bond_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "bond delta")?;
    match op {
        "add" => Ok(DeltaInput::BondAdd(parse_bond_entry(payload)?)),
        "remove" => Ok(DeltaInput::BondRemove(BondRef::from_edn(payload)?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "bond :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["bond delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "bond :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::BondModify(
                BondRef::from_edn(&v[0])?,
                BondUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!("unknown bond delta op :{o}"))),
    }
}

fn parse_stereo_kind(edn: &Edn<'_>) -> Result<StereoKind, DeError> {
    match edn {
        Edn::Keyword(k) => stereo_kind_from_name(k.name()),
        other => Err(DeError::TypeMismatch {
            expected: "stereo kind keyword",
            got: other.kind(),
            path: Vec::new(),
        }),
    }
}

fn parse_delta_dative_bond_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "dative-bond delta")?;
    match op {
        "add" => Ok(DeltaInput::DativeBondAdd(parse_dative_bond_entry(payload)?)),
        "remove" => Ok(DeltaInput::DativeBondRemove(DativeBondRef::from_edn(
            payload,
        )?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "dative-bond :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["dative-bond delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "dative-bond :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::DativeBondModify(
                DativeBondRef::from_edn(&v[0])?,
                DativeBondUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!(
            "unknown dative-bond delta op :{o}"
        ))),
    }
}

fn parse_delta_aromatic_system_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "aromatic-system delta")?;
    match op {
        "add" => Ok(DeltaInput::AromaticSystemAdd(parse_aromatic_system_entry(
            payload,
        )?)),
        "remove" => Ok(DeltaInput::AromaticSystemRemove(
            AromaticSystemRef::from_edn(payload)?,
        )),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "aromatic-system :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["aromatic-system delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "aromatic-system :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::AromaticSystemModify(
                AromaticSystemRef::from_edn(&v[0])?,
                AromaticSystemUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!(
            "unknown aromatic-system delta op :{o}"
        ))),
    }
}

fn parse_delta_multicenter_bond_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "multicenter-bond delta")?;
    match op {
        "add" => Ok(DeltaInput::MulticenterBondAdd(
            parse_multicenter_bond_entry(payload)?,
        )),
        "remove" => Ok(DeltaInput::MulticenterBondRemove(
            MulticenterBondRef::from_edn(payload)?,
        )),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "multicenter-bond :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["multicenter-bond delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "multicenter-bond :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::MulticenterBondModify(
                MulticenterBondRef::from_edn(&v[0])?,
                MulticenterBondUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!(
            "unknown multicenter-bond delta op :{o}"
        ))),
    }
}

fn parse_delta_noncovalent_bond_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "noncovalent-bond delta")?;
    match op {
        "add" => Ok(DeltaInput::NoncovalentBondAdd(
            parse_noncovalent_bond_entry(payload)?,
        )),
        "remove" => Ok(DeltaInput::NoncovalentBondRemove(
            NoncovalentBondRef::from_edn(payload)?,
        )),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "noncovalent-bond :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["noncovalent-bond delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "noncovalent-bond :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::NoncovalentBondModify(
                NoncovalentBondRef::from_edn(&v[0])?,
                NoncovalentBondUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!(
            "unknown noncovalent-bond delta op :{o}"
        ))),
    }
}

fn parse_delta_stereo_atom_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "stereo-atom delta")?;
    match op {
        "add" => Ok(DeltaInput::StereoAtomAdd(parse_stereo_atom_entry(payload)?)),
        "remove" => Ok(DeltaInput::StereoAtomRemove(StereoAtomRef::from_edn(
            payload,
        )?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "stereo-atom :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["stereo-atom delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "stereo-atom :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::StereoAtomModify(
                StereoAtomRef::from_edn(&v[0])?,
                StereoAtomUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        "swap" => {
            let (r, kind) = parse_stereo_transform(payload, "stereo-atom :swap")?;
            Ok(DeltaInput::StereoAtomSwap(
                StereoAtomRef::from_edn(r)?,
                kind,
            ))
        }
        "mirror" => {
            let (r, kind) = parse_stereo_transform(payload, "stereo-atom :mirror")?;
            Ok(DeltaInput::StereoAtomMirror(
                StereoAtomRef::from_edn(r)?,
                kind,
            ))
        }
        "apply" => {
            let (r, kind, permutation) = parse_stereo_apply(payload, "stereo-atom :apply")?;
            Ok(DeltaInput::StereoAtomApply(
                StereoAtomRef::from_edn(r)?,
                kind,
                permutation,
            ))
        }
        o => Err(DeError::Custom(format!(
            "unknown stereo-atom delta op :{o}"
        ))),
    }
}

fn parse_delta_stereo_bond_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "stereo-bond delta")?;
    match op {
        "add" => Ok(DeltaInput::StereoBondAdd(parse_stereo_bond_entry(payload)?)),
        "remove" => Ok(DeltaInput::StereoBondRemove(StereoBondRef::from_edn(
            payload,
        )?)),
        "modify" => {
            let Edn::Vector(v) = payload else {
                return Err(DeError::TypeMismatch {
                    expected: "stereo-bond :modify [ref dsl]",
                    got: payload.kind(),
                    path: vec!["stereo-bond delta".into()],
                });
            };
            if v.len() != 2 {
                return Err(DeError::Custom(format!(
                    "stereo-bond :modify expects [ref dsl], got {} elements",
                    v.len()
                )));
            }
            Ok(DeltaInput::StereoBondModify(
                StereoBondRef::from_edn(&v[0])?,
                StereoBondUpdateDsl::from_edn(&v[1])?.0,
            ))
        }
        "swap" => {
            let (r, kind) = parse_stereo_transform(payload, "stereo-bond :swap")?;
            Ok(DeltaInput::StereoBondSwap(
                StereoBondRef::from_edn(r)?,
                kind,
            ))
        }
        "mirror" => {
            let (r, kind) = parse_stereo_transform(payload, "stereo-bond :mirror")?;
            Ok(DeltaInput::StereoBondMirror(
                StereoBondRef::from_edn(r)?,
                kind,
            ))
        }
        "apply" => {
            let (r, kind, permutation) = parse_stereo_apply(payload, "stereo-bond :apply")?;
            Ok(DeltaInput::StereoBondApply(
                StereoBondRef::from_edn(r)?,
                kind,
                permutation,
            ))
        }
        o => Err(DeError::Custom(format!(
            "unknown stereo-bond delta op :{o}"
        ))),
    }
}

/// Extract the `[<ref> <kind>]` payload shared by the `:swap` / `:mirror` stereo verbs; the
/// ref element is returned unresolved for the caller's ref type.
fn parse_stereo_transform<'a>(
    payload: &'a Edn<'a>,
    ctx: &'static str,
) -> Result<(&'a Edn<'a>, StereoKind), DeError> {
    let Edn::Vector(v) = payload else {
        return Err(DeError::TypeMismatch {
            expected: ctx,
            got: payload.kind(),
            path: Vec::new(),
        });
    };
    if v.len() != 2 {
        return Err(DeError::Custom(format!(
            "{ctx} expects [ref kind], got {} elements",
            v.len()
        )));
    }
    Ok((&v[0], parse_stereo_kind(&v[1])?))
}

/// Extract the `[<ref> <kind> <cycles>]` payload for the `:apply` stereo verb; the permutation
/// degree comes from the explicit kind, not the referenced entity.
fn parse_stereo_apply<'a>(
    payload: &'a Edn<'a>,
    ctx: &'static str,
) -> Result<(&'a Edn<'a>, StereoKind, Permutation), DeError> {
    let Edn::Vector(v) = payload else {
        return Err(DeError::TypeMismatch {
            expected: ctx,
            got: payload.kind(),
            path: Vec::new(),
        });
    };
    if v.len() != 3 {
        return Err(DeError::Custom(format!(
            "{ctx} expects [ref kind cycles], got {} elements",
            v.len()
        )));
    }
    let kind = parse_stereo_kind(&v[1])?;
    let Edn::Str(s) = &v[2] else {
        return Err(DeError::TypeMismatch {
            expected: "permutation cycle string",
            got: v[2].kind(),
            path: Vec::new(),
        });
    };
    let permutation =
        parse_permutation(s.as_ref(), kind.degree()).map_err(|e| DeError::subgrammar(ctx, e))?;
    Ok((&v[0], kind, permutation))
}

fn parse_delta_constraint_input(edn: &Edn<'_>) -> Result<DeltaInput, DeError> {
    let (op, payload) = parse_single_key_map(edn, "constraint delta")?;
    match op {
        "add" => Ok(DeltaInput::ConstraintAdd(ConstraintDsl::from_edn(payload)?)),
        "remove" => Ok(DeltaInput::ConstraintRemove(ConstraintDsl::from_edn(
            payload,
        )?)),
        o => Err(DeError::Custom(format!("unknown constraint delta op :{o}"))),
    }
}

/// Render a reaction to its EDN map: `:lhs` (the molecule via the molecule renderer), `:deltas`, and
/// `:atom-aliases` (reaction-level, when present). Aliases render last, as in the molecule surface.
fn render_reaction_edn(ast: &ReactionAst, meta: &ReactionMetadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(3);
    map.insert(
        Edn::keyword("lhs"),
        render_molecule_edn(&ast.lhs, meta.lhs()),
    );
    map.insert(
        Edn::keyword("deltas"),
        Edn::Vector(render_deltas(&ast.deltas, meta).into()),
    );
    if meta.has_atom_aliases() {
        let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(meta.atom_aliases_len() * 2);
        for (name, dsl) in meta.iter_atom_aliases() {
            pairs.push(Edn::Keyword(EdnKeyword::owned(name.to_string())));
            pairs.push(dsl.to_edn());
        }
        map.insert(Edn::keyword("atom-aliases"), Edn::Vector(pairs.into()));
    }
    Edn::Map(map)
}

/// Render `deltas` to the `:deltas` vector entries, resolving refs against the reaction `meta`.
/// Consecutive `ModifyField` / `ModifyConstraint` for one atom coalesce into a single
/// `{:atom {:modify [<ref> <partial>]}}`.
fn render_deltas(deltas: &Deltas, meta: &ReactionMetadata) -> Vec<Edn<'static>> {
    let deltas = deltas.as_slice();
    let mut out = Vec::new();
    // Built once on the first constraint delta (constraint refs resolve against lhs ∪ created).
    let mut combined_metadata: Option<MoleculeMetadata> = None;
    let mut i = 0;
    while i < deltas.len() {
        match &deltas[i] {
            Delta::Atom(AtomDelta::Add { id, ast }) => {
                out.push(single_key_map(
                    "atom",
                    single_key_map("add", render_atom_entry(*id, ast, meta)),
                ));
                i += 1;
            }
            Delta::Atom(AtomDelta::Remove { id, .. }) => {
                out.push(single_key_map(
                    "atom",
                    single_key_map("remove", render_atom_ref(*id, meta)),
                ));
                i += 1;
            }
            Delta::Atom(
                AtomDelta::ModifyField { id, .. } | AtomDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = AtomUpdate::default();
                while let Some(Delta::Atom(delta)) = deltas.get(i) {
                    match delta {
                        AtomDelta::ModifyField { id: j, change } if *j == id => match change {
                            AtomFieldChange::Element { new, .. } => {
                                update.element = Some(new.clone())
                            }
                            AtomFieldChange::IsotopeMass { new, .. } => {
                                update.isotope_mass = Some(new.clone())
                            }
                            AtomFieldChange::Charge { new, .. } => {
                                update.charge = Some(new.clone())
                            }
                            AtomFieldChange::ImplicitHydrogens { new, .. } => {
                                update.implicit_hydrogens = Some(new.clone())
                            }
                            AtomFieldChange::LonePairs { new, .. } => {
                                update.lone_pairs = Some(new.clone())
                            }
                            AtomFieldChange::Spin { new, .. } => {
                                update.spin.unpaired = Some(new.unpaired.clone());
                                update.spin.multiplicity = Some(new.multiplicity.clone());
                            }
                        },
                        AtomDelta::ModifyConstraint { id: j, old, new } if *j == id => match new {
                            Some(c) => {
                                update.constraints.set(c.clone());
                            }
                            None => {
                                if let Some(old) = old {
                                    update.constraints.set(old.as_undetermined());
                                }
                            }
                        },
                        _ => break,
                    }
                    i += 1;
                }
                let payload = Edn::Vector(
                    vec![render_atom_ref(id, meta), AtomUpdateDsl(update).to_edn()].into(),
                );
                out.push(single_key_map("atom", single_key_map("modify", payload)));
            }
            Delta::Bond(BondDelta::Add { id, atoms, ast }) => {
                out.push(single_key_map(
                    "bond",
                    single_key_map("add", render_bond_entry(*id, *atoms, ast, meta)),
                ));
                i += 1;
            }
            Delta::Bond(BondDelta::Remove { id, .. }) => {
                out.push(single_key_map(
                    "bond",
                    single_key_map("remove", render_bond_ref(*id, meta)),
                ));
                i += 1;
            }
            Delta::Bond(
                BondDelta::ModifyField { id, .. } | BondDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = BondUpdate::default();
                while let Some(Delta::Bond(delta)) = deltas.get(i) {
                    match delta {
                        BondDelta::ModifyField { id: j, change } if *j == id => match change {
                            BondFieldChange::Order { new, .. } => update.order = Some(new.clone()),
                            BondFieldChange::Charge { new, .. } => {
                                update.charge = Some(new.clone())
                            }
                            BondFieldChange::Spin { new, .. } => {
                                update.spin.unpaired = Some(new.unpaired.clone());
                                update.spin.multiplicity = Some(new.multiplicity.clone());
                            }
                        },
                        BondDelta::ModifyConstraint { id: j, old, new } if *j == id => match new {
                            Some(c) => {
                                update.constraints.set(c.clone());
                            }
                            None => {
                                if let Some(old) = old {
                                    update.constraints.set(old.as_undetermined());
                                }
                            }
                        },
                        _ => break,
                    }
                    i += 1;
                }
                let payload = Edn::Vector(
                    vec![render_bond_ref(id, meta), BondUpdateDsl(update).to_edn()].into(),
                );
                out.push(single_key_map("bond", single_key_map("modify", payload)));
            }
            Delta::Constraint(delta) => {
                let (op, constraint) = match delta {
                    ConstraintDelta::Add(c) => ("add", c),
                    ConstraintDelta::Remove(c) => ("remove", c),
                };
                let combined = combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let dsl = ConstraintDsl::from_ast(constraint, combined)
                    .expect("ConstraintDsl::from_ast is infallible for a well-formed AST");
                out.push(single_key_map(
                    "constraint",
                    single_key_map(op, dsl.to_edn()),
                ));
                i += 1;
            }
            Delta::DativeBond(DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ast,
            }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "dative-bond",
                    single_key_map(
                        "add",
                        render_dative_entry(
                            *id,
                            donors.iter().copied(),
                            *acceptor,
                            DativeBondDsl::from_ref(ast).to_edn(),
                            combined,
                        ),
                    ),
                ));
                i += 1;
            }
            Delta::DativeBond(DativeBondDelta::Remove { id, .. }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "dative-bond",
                    single_key_map("remove", DativeBondRef::denote(*id, combined).to_edn()),
                ));
                i += 1;
            }
            Delta::DativeBond(
                DativeBondDelta::ModifyField { id, .. }
                | DativeBondDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = DativeBondUpdate::default();
                while let Some(Delta::DativeBond(delta)) = deltas.get(i) {
                    match delta {
                        DativeBondDelta::ModifyField { id: j, change } if *j == id => {
                            match change {
                                DativeBondFieldChange::Order { new, .. } => {
                                    update.order = Some(new.clone())
                                }
                            }
                        }
                        DativeBondDelta::ModifyConstraint { id: j, old, new } if *j == id => {
                            match new {
                                Some(c) => {
                                    update.constraints.set(c.clone());
                                }
                                None => {
                                    if let Some(old) = old {
                                        update.constraints.set(old.as_undetermined());
                                    }
                                }
                            }
                        }
                        _ => break,
                    }
                    i += 1;
                }
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        DativeBondRef::denote(id, combined).to_edn(),
                        DativeBondUpdateDsl(update).to_edn(),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "dative-bond",
                    single_key_map("modify", payload),
                ));
            }
            Delta::AromaticSystem(AromaticSystemDelta::Add { id, atoms, ast }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "aromatic-system",
                    single_key_map(
                        "add",
                        render_aromatic_entry(
                            *id,
                            atoms.iter().copied(),
                            AromaticSystemDsl::from_ref(ast).to_edn(),
                            combined,
                        ),
                    ),
                ));
                i += 1;
            }
            Delta::AromaticSystem(AromaticSystemDelta::Remove { id, .. }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "aromatic-system",
                    single_key_map("remove", AromaticSystemRef::denote(*id, combined).to_edn()),
                ));
                i += 1;
            }
            Delta::AromaticSystem(
                AromaticSystemDelta::ModifyField { id, .. }
                | AromaticSystemDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = AromaticSystemUpdate::default();
                while let Some(Delta::AromaticSystem(delta)) = deltas.get(i) {
                    match delta {
                        AromaticSystemDelta::ModifyField { id: j, change } if *j == id => {
                            match change {
                                AromaticSystemFieldChange::Electrons { new, .. } => {
                                    update.electrons = Some(new.clone())
                                }
                                AromaticSystemFieldChange::Charge { new, .. } => {
                                    update.charge = Some(new.clone())
                                }
                                AromaticSystemFieldChange::Spin { new, .. } => {
                                    update.spin.unpaired = Some(new.unpaired.clone());
                                    update.spin.multiplicity = Some(new.multiplicity.clone());
                                }
                            }
                        }
                        AromaticSystemDelta::ModifyConstraint { id: j, old, new } if *j == id => {
                            match new {
                                Some(c) => {
                                    update.constraints.set(c.clone());
                                }
                                None => {
                                    if let Some(old) = old {
                                        update.constraints.set(old.as_undetermined());
                                    }
                                }
                            }
                        }
                        _ => break,
                    }
                    i += 1;
                }
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        AromaticSystemRef::denote(id, combined).to_edn(),
                        AromaticSystemUpdateDsl(update).to_edn(),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "aromatic-system",
                    single_key_map("modify", payload),
                ));
            }
            Delta::MulticenterBond(MulticenterBondDelta::Add { id, atoms, ast }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "multicenter-bond",
                    single_key_map(
                        "add",
                        render_multicenter_entry(
                            *id,
                            atoms.iter().copied(),
                            MulticenterBondDsl::from_ref(ast).to_edn(),
                            combined,
                        ),
                    ),
                ));
                i += 1;
            }
            Delta::MulticenterBond(MulticenterBondDelta::Remove { id, .. }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "multicenter-bond",
                    single_key_map("remove", MulticenterBondRef::denote(*id, combined).to_edn()),
                ));
                i += 1;
            }
            Delta::MulticenterBond(
                MulticenterBondDelta::ModifyField { id, .. }
                | MulticenterBondDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = MulticenterBondUpdate::default();
                while let Some(Delta::MulticenterBond(delta)) = deltas.get(i) {
                    match delta {
                        MulticenterBondDelta::ModifyField { id: j, change } if *j == id => {
                            match change {
                                MulticenterBondFieldChange::Electrons { new, .. } => {
                                    update.electrons = Some(new.clone())
                                }
                                MulticenterBondFieldChange::Charge { new, .. } => {
                                    update.charge = Some(new.clone())
                                }
                                MulticenterBondFieldChange::Spin { new, .. } => {
                                    update.spin.unpaired = Some(new.unpaired.clone());
                                    update.spin.multiplicity = Some(new.multiplicity.clone());
                                }
                            }
                        }
                        MulticenterBondDelta::ModifyConstraint { id: j, old, new } if *j == id => {
                            match new {
                                Some(c) => {
                                    update.constraints.set(c.clone());
                                }
                                None => {
                                    if let Some(old) = old {
                                        update.constraints.set(old.as_undetermined());
                                    }
                                }
                            }
                        }
                        _ => break,
                    }
                    i += 1;
                }
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        MulticenterBondRef::denote(id, combined).to_edn(),
                        MulticenterBondUpdateDsl(update).to_edn(),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "multicenter-bond",
                    single_key_map("modify", payload),
                ));
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Add { id, atoms, ast }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "noncovalent-bond",
                    single_key_map(
                        "add",
                        render_noncovalent_entry(
                            *id,
                            *atoms,
                            NoncovalentBondDsl::from_ref(ast).to_edn(),
                            combined,
                        ),
                    ),
                ));
                i += 1;
            }
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, .. }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "noncovalent-bond",
                    single_key_map("remove", NoncovalentBondRef::denote(*id, combined).to_edn()),
                ));
                i += 1;
            }
            Delta::NoncovalentBond(
                NoncovalentBondDelta::ModifyField { id, .. }
                | NoncovalentBondDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = NoncovalentBondUpdate::default();
                while let Some(Delta::NoncovalentBond(delta)) = deltas.get(i) {
                    match delta {
                        NoncovalentBondDelta::ModifyField { id: j, change } if *j == id => {
                            match change {
                                NoncovalentBondFieldChange::Kind { new, .. } => {
                                    update.kind = Some(new.clone())
                                }
                            }
                        }
                        NoncovalentBondDelta::ModifyConstraint { id: j, old, new } if *j == id => {
                            match new {
                                Some(c) => {
                                    update.constraints.set(c.clone());
                                }
                                None => {
                                    if let Some(old) = old {
                                        update.constraints.set(old.as_undetermined());
                                    }
                                }
                            }
                        }
                        _ => break,
                    }
                    i += 1;
                }
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        NoncovalentBondRef::denote(id, combined).to_edn(),
                        NoncovalentBondUpdateDsl(update).to_edn(),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "noncovalent-bond",
                    single_key_map("modify", payload),
                ));
            }
            Delta::StereoAtom(StereoAtomDelta::Add {
                id,
                site,
                ligands,
                ast,
            }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "stereo-atom",
                    single_key_map(
                        "add",
                        render_stereo_atom_entry(
                            *id,
                            *site,
                            ligands
                                .iter()
                                .map(|l| render_stereo_ligand(*l, combined))
                                .collect(),
                            StereoAtomDsl::from_ref(ast).to_edn(),
                            combined,
                        ),
                    ),
                ));
                i += 1;
            }
            Delta::StereoAtom(StereoAtomDelta::Remove { id, .. }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "stereo-atom",
                    single_key_map("remove", StereoAtomRef::denote(*id, combined).to_edn()),
                ));
                i += 1;
            }
            Delta::StereoAtom(
                StereoAtomDelta::ModifyField { id, .. }
                | StereoAtomDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = StereoAtomUpdate::default();
                let mut modify_kind: Option<StereoKind> = None;
                while let Some(Delta::StereoAtom(delta)) = deltas.get(i) {
                    match delta {
                        StereoAtomDelta::ModifyField { id: j, change } if *j == id => {
                            let clears_configuration = match change {
                                StereoAtomFieldChange::Configuration { new, .. } => {
                                    update.configuration = match new {
                                        StereoConfigurationAst::Undetermined => {
                                            StereoConfigurationUpdate::Undetermined
                                        }
                                        StereoConfigurationAst::Kinded(kind, coset) => {
                                            StereoConfigurationUpdate::Kinded {
                                                kind: *kind,
                                                coset: Some(coset.clone()),
                                            }
                                        }
                                    };
                                    matches!(new, StereoConfigurationAst::Undetermined)
                                }
                            };
                            i += 1;
                            if clears_configuration {
                                break;
                            }
                            continue;
                        }
                        StereoAtomDelta::ModifyConstraint {
                            id: j,
                            kind,
                            old,
                            new,
                        } if *j == id => {
                            modify_kind = modify_kind.or(*kind);
                            match new {
                                Some(c) => {
                                    update.constraints.set(c.clone());
                                }
                                None => {
                                    if let Some(old) = old {
                                        update.constraints.set(old.as_undetermined());
                                    }
                                }
                            }
                        }
                        _ => break,
                    }
                    i += 1;
                }
                if matches!(update.configuration, StereoConfigurationUpdate::Unchanged) {
                    if let Some(kind) = modify_kind {
                        update.configuration =
                            StereoConfigurationUpdate::Kinded { kind, coset: None };
                    }
                }
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoAtomRef::denote(id, combined).to_edn(),
                        StereoAtomUpdateDsl(update).to_edn(),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-atom",
                    single_key_map("modify", payload),
                ));
            }
            Delta::StereoAtom(StereoAtomDelta::Swap { id, kind }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoAtomRef::denote(*id, combined).to_edn(),
                        render_edn_stereo_kind(*kind),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-atom",
                    single_key_map("swap", payload),
                ));
                i += 1;
            }
            Delta::StereoAtom(StereoAtomDelta::Mirror { id, kind }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoAtomRef::denote(*id, combined).to_edn(),
                        render_edn_stereo_kind(*kind),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-atom",
                    single_key_map("mirror", payload),
                ));
                i += 1;
            }
            Delta::StereoAtom(StereoAtomDelta::Apply {
                id,
                kind,
                permutation,
            }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoAtomRef::denote(*id, combined).to_edn(),
                        render_edn_stereo_kind(*kind),
                        Edn::Str(Cow::Owned(permutation.to_string())),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-atom",
                    single_key_map("apply", payload),
                ));
                i += 1;
            }
            Delta::StereoBond(StereoBondDelta::Add {
                id,
                site,
                ligands,
                ast,
            }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "stereo-bond",
                    single_key_map(
                        "add",
                        render_stereo_bond_entry(
                            *id,
                            *site,
                            ligands
                                .iter()
                                .map(|l| render_stereo_ligand(*l, combined))
                                .collect(),
                            StereoBondDsl::from_ref(ast).to_edn(),
                            combined,
                        ),
                    ),
                ));
                i += 1;
            }
            Delta::StereoBond(StereoBondDelta::Remove { id, .. }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                out.push(single_key_map(
                    "stereo-bond",
                    single_key_map("remove", StereoBondRef::denote(*id, combined).to_edn()),
                ));
                i += 1;
            }
            Delta::StereoBond(
                StereoBondDelta::ModifyField { id, .. }
                | StereoBondDelta::ModifyConstraint { id, .. },
            ) => {
                let id = *id;
                let mut update = StereoBondUpdate::default();
                let mut modify_kind: Option<StereoKind> = None;
                while let Some(Delta::StereoBond(delta)) = deltas.get(i) {
                    match delta {
                        StereoBondDelta::ModifyField { id: j, change } if *j == id => {
                            let clears_configuration = match change {
                                StereoBondFieldChange::Configuration { new, .. } => {
                                    update.configuration = match new {
                                        StereoConfigurationAst::Undetermined => {
                                            StereoConfigurationUpdate::Undetermined
                                        }
                                        StereoConfigurationAst::Kinded(kind, coset) => {
                                            StereoConfigurationUpdate::Kinded {
                                                kind: *kind,
                                                coset: Some(coset.clone()),
                                            }
                                        }
                                    };
                                    matches!(new, StereoConfigurationAst::Undetermined)
                                }
                            };
                            i += 1;
                            if clears_configuration {
                                break;
                            }
                            continue;
                        }
                        StereoBondDelta::ModifyConstraint {
                            id: j,
                            kind,
                            old,
                            new,
                        } if *j == id => {
                            modify_kind = modify_kind.or(*kind);
                            match new {
                                Some(c) => {
                                    update.constraints.set(c.clone());
                                }
                                None => {
                                    if let Some(old) = old {
                                        update.constraints.set(old.as_undetermined());
                                    }
                                }
                            }
                        }
                        _ => break,
                    }
                    i += 1;
                }
                if matches!(update.configuration, StereoConfigurationUpdate::Unchanged) {
                    if let Some(kind) = modify_kind {
                        update.configuration =
                            StereoConfigurationUpdate::Kinded { kind, coset: None };
                    }
                }
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoBondRef::denote(id, combined).to_edn(),
                        StereoBondUpdateDsl(update).to_edn(),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-bond",
                    single_key_map("modify", payload),
                ));
            }
            Delta::StereoBond(StereoBondDelta::Swap { id, kind }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoBondRef::denote(*id, combined).to_edn(),
                        render_edn_stereo_kind(*kind),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-bond",
                    single_key_map("swap", payload),
                ));
                i += 1;
            }
            Delta::StereoBond(StereoBondDelta::Mirror { id, kind }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoBondRef::denote(*id, combined).to_edn(),
                        render_edn_stereo_kind(*kind),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-bond",
                    single_key_map("mirror", payload),
                ));
                i += 1;
            }
            Delta::StereoBond(StereoBondDelta::Apply {
                id,
                kind,
                permutation,
            }) => {
                let combined = &*combined_metadata.get_or_insert_with(|| meta.combined_metadata());
                let payload = Edn::Vector(
                    vec![
                        StereoBondRef::denote(*id, combined).to_edn(),
                        render_edn_stereo_kind(*kind),
                        Edn::Str(Cow::Owned(permutation.to_string())),
                    ]
                    .into(),
                );
                out.push(single_key_map(
                    "stereo-bond",
                    single_key_map("apply", payload),
                ));
                i += 1;
            }
        }
    }
    out
}

/// A delta ref (`:remove` / `:modify`) names an existing lhs entity — resolved against the lhs frame.
fn render_atom_ref(id: AtomId, meta: &ReactionMetadata) -> Edn<'static> {
    match meta.lhs().keyword(Entity::Atom(id)) {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A created atom (`:add`) — its keyword and aliases live in the reaction frame (the alias namespace
/// is the lhs ∪ reaction union). Renders `<atom-dsl>` or `[<keyword> <atom-dsl>]`.
fn render_atom_entry(id: AtomId, atom: &AtomAst, meta: &ReactionMetadata) -> Edn<'static> {
    let dsl = AtomDsl::from_ref(atom);
    let spec = match meta
        .atom_alias_for(dsl)
        .or_else(|| meta.lhs().atom_alias_for(dsl))
    {
        Some(alias) => Edn::Keyword(EdnKeyword::owned(alias.to_string())),
        None => dsl.to_edn(),
    };
    match meta.atom_keyword(id) {
        Some(name) => {
            Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(name.to_string())), spec].into())
        }
        None => spec,
    }
}

/// An atom named as a bond endpoint — resolved against the union namespace (lhs ∪ created), since a
/// bond may attach to a same-reaction atom. Unlike a delta target ref, which is lhs-only.
fn render_atom_endpoint(id: AtomId, meta: &ReactionMetadata) -> Edn<'static> {
    match meta
        .atom_keyword(id)
        .or_else(|| meta.lhs().keyword(Entity::Atom(id)))
    {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A bond delta target (`:remove` / `:modify`) names an existing lhs bond — resolved lhs-frame only.
fn render_bond_ref(id: BondId, meta: &ReactionMetadata) -> Edn<'static> {
    match meta.lhs().keyword(Entity::Bond(id)) {
        Some(name) => Edn::Keyword(EdnKeyword::owned(name.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A created bond (`:add`). Renders `[<a> <b> <bond-dsl>]`, or
/// `{:id <keyword> :atoms [<a> <b>] :type <bond-dsl>}` when the bond carries a keyword.
/// Endpoints resolve against the union namespace.
fn render_bond_entry(
    id: BondId,
    atoms: [AtomId; 2],
    ast: &BondAst,
    meta: &ReactionMetadata,
) -> Edn<'static> {
    let bond_edn = BondDsl::from_ref(ast).to_edn();
    let first = render_atom_endpoint(atoms[0], meta);
    let second = render_atom_endpoint(atoms[1], meta);
    match meta.bond_keyword(id) {
        Some(name) => {
            let mut m = EdnMap::with_capacity(3);
            m.insert(
                Edn::keyword("id"),
                Edn::Keyword(EdnKeyword::owned(name.to_string())),
            );
            m.insert(
                Edn::keyword("atoms"),
                Edn::Vector(vec![first, second].into()),
            );
            m.insert(Edn::keyword("type"), bond_edn);
            Edn::Map(m)
        }
        None => Edn::Vector(vec![first, second, bond_edn].into()),
    }
}

#[cfg(test)]
mod tests {
    use std::fs;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::atom::ElementAst;
    use crate::ast::boolean::BooleanAst;
    use crate::ast::constraint::{
        AromaticSystemConstraintAst, AromaticSystemConstraintsAst, AtomConstraintAst,
        BondConstraintAst, Constraint, DativeBondConstraintAst, DativeBondConstraintsAst,
        MoleculeConstraint, MulticenterBondConstraintAst, MulticenterBondConstraintsAst,
        NoncovalentBondConstraintAst, NoncovalentBondConstraintsAst, RingScope,
        StereoAtomConstraintAst, StereoAtomConstraintsAst, StereoBondConstraintAst,
        StereoBondConstraintsAst, StereogenicityAst,
    };
    use crate::ast::delta::{ConstraintDelta, Deltas};
    use crate::ast::edit::{
        AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, MulticenterBondFieldChange,
        NoncovalentBondFieldChange,
    };
    use crate::ast::electrons::ElectronCountsAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::noncovalent::{NoncovalentBondKind, NoncovalentBondKindAst};
    use crate::ast::spin::{SpinStateAst, SpinStateUpdate};
    use crate::ast::stereo::{StereoCoset, Stereogenicity};
    use crate::ast::value::ValueAst;
    use crate::dsl::bond::BondDsl;
    use crate::dsl::constraint::MoleculeConstraintDsl;
    use crate::dsl::molecule::AtomSpecInput;
    use crate::dsl::refs::{
        AromaticSystemRef, AtomRef, BondRef, DativeBondRef, MulticenterBondRef, NoncovalentBondRef,
    };
    use crate::mol_dsl;

    #[rstest]
    #[case::sn2(ReactionAst::new(
        mol_dsl!(r##"{:atoms ["C" "Br"] :bonds [[0 1 "1"]]}"##),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Add { id: AtomId(2), ast: AtomAst::from_element(Element::O) }),
            Delta::Bond(BondDelta::Add {
                id: BondId(1),
                atoms: [AtomId(0), AtomId(2)],
                ast: BondAst::from_order(1),
            }),
            Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(1),
                change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) },
            }),
            Delta::Constraint(ConstraintDelta::Add(Constraint::And(vec![]))),
        ]),
    ))]
    #[case::dative(r##"{:lhs {:atoms ["C" "N"]} :deltas [{:dative-bond {:add {:donors [0] :acceptor 1 :type "1#R"}}}]}"##.parse().unwrap())]
    #[case::aromatic(r##"{:lhs {:atoms ["C" "C"]} :deltas [{:aromatic-system {:add {:atoms [0 1] :type "*#e2"}}}]}"##.parse().unwrap())]
    #[case::multicenter(r##"{:lhs {:atoms ["B" "H" "B"]} :deltas [{:multicenter-bond {:add {:atoms [0 1 2] :type "[1,0,1]#e2"}}}]}"##.parse().unwrap())]
    #[case::noncovalent(r##"{:lhs {:atoms ["N" "H"]} :deltas [{:noncovalent-bond {:add {:atoms [0 1] :type "Hbd"}}}]}"##.parse().unwrap())]
    #[case::stereo_atom(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]} :deltas [{:stereo-atom {:add {:site 0 :ligands [1 2 3 4] :type "Th1"}}}]}"##.parse().unwrap())]
    #[case::stereo_bond(r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]} :deltas [{:stereo-bond {:add {:site 1 :ligands [0 3] :type "Ct1"}}}]}"##.parse().unwrap())]
    fn test_reaction_dsl_from_ast_roundtrip(#[case] reaction: ReactionAst) {
        let cfg = ReactionDefaults::ground();
        let dsl = ReactionDsl::from_parts(reaction, ReactionMetadata::default());
        let lowered = ReactionDsl::from_ast(&dsl.clone().into_ast(&cfg), &cfg);
        assert_eq!(lowered.ast(), dsl.ast());
    }

    #[rstest]
    #[case::add_bare(
        r##"{:atom {:add "C#h3"}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            keyword: None,
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                a.implicit_hydrogens = ValueAst::Lit(3);
                a
            }))),
        })
    )]
    #[case::add_keyword(
        r##"{:atom {:add [:nu "O#h1"]}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            keyword: Some("nu".into()),
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::O));
                a.implicit_hydrogens = ValueAst::Lit(1);
                a
            }))),
        })
    )]
    #[case::add_alias(
        "{:atom {:add :foo}}",
        DeltaInput::AtomAdd(AtomEntryInput {
            keyword: None,
            spec: AtomSpecInput::Alias("foo".into()),
        })
    )]
    #[case::remove_keyword("{:atom {:remove :br}}", DeltaInput::AtomRemove(AtomRef::Keyword("br".into())))]
    #[case::remove_index("{:atom {:remove 1}}", DeltaInput::AtomRemove(AtomRef::Index(1)))]
    #[case::modify(
        r##"{:atom {:modify [:br "#c-1"]}}"##,
        DeltaInput::AtomModify(
            AtomRef::Keyword("br".into()),
            AtomUpdate { charge: Some(ValueAst::Lit(-1)), ..Default::default() },
        )
    )]
    fn test_parse_delta_input_atom(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add_bare(
        r##"{:atom {:add "C#h3"}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            keyword: None,
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                a.implicit_hydrogens = ValueAst::Lit(3);
                a
            }))),
        })
    )]
    #[case::add_keyword(
        r##"{:atom {:add [:nu "O#h1"]}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            keyword: Some("nu".into()),
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::O));
                a.implicit_hydrogens = ValueAst::Lit(1);
                a
            }))),
        })
    )]
    #[case::add_alias(
        "{:atom {:add :foo}}",
        DeltaInput::AtomAdd(AtomEntryInput {
            keyword: None,
            spec: AtomSpecInput::Alias("foo".into()),
        })
    )]
    #[case::remove_keyword("{:atom {:remove :br}}", DeltaInput::AtomRemove(AtomRef::Keyword("br".into())))]
    #[case::remove_index("{:atom {:remove 1}}", DeltaInput::AtomRemove(AtomRef::Index(1)))]
    #[case::modify(
        r##"{:atom {:modify [:br "#c-1"]}}"##,
        DeltaInput::AtomModify(
            AtomRef::Keyword("br".into()),
            AtomUpdate { charge: Some(ValueAst::Lit(-1)), ..Default::default() },
        )
    )]
    fn test_read_delta_input_atom(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add_vec(
        r##"{:bond {:add [0 1 "1"]}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            keyword: None,
            first: AtomRef::Index(0),
            second: AtomRef::Index(1),
            bond: BondDsl(BondAst::from_order(1)),
        })
    )]
    #[case::add_map_keyword(
        r##"{:bond {:add {:id :b1 :atoms [:c :nu] :type "2"}}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            keyword: Some("b1".into()),
            first: AtomRef::Keyword("c".into()),
            second: AtomRef::Keyword("nu".into()),
            bond: BondDsl(BondAst::from_order(2)),
        })
    )]
    #[case::remove_keyword("{:bond {:remove :b1}}", DeltaInput::BondRemove(BondRef::Keyword("b1".into())))]
    #[case::remove_index("{:bond {:remove 0}}", DeltaInput::BondRemove(BondRef::Index(0)))]
    #[case::modify(
        r##"{:bond {:modify [:b1 "2"]}}"##,
        DeltaInput::BondModify(
            BondRef::Keyword("b1".into()),
            BondUpdate { order: Some(ValueAst::Lit(2)), ..Default::default() },
        )
    )]
    #[case::modify_undetermined(
        r##"{:bond {:modify [:b1 "*#c*#u*#s*"]}}"##,
        DeltaInput::BondModify(
            BondRef::Keyword("b1".into()),
            BondUpdate {
                order: Some(ValueAst::Undetermined),
                charge: Some(ValueAst::Undetermined),
                spin: SpinStateUpdate {
                    unpaired: Some(ValueAst::Undetermined),
                    multiplicity: Some(ValueAst::Undetermined),
                },
                ..Default::default()
            },
        )
    )]
    fn test_parse_delta_input_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add_vec(
        r##"{:bond {:add [0 1 "1"]}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            keyword: None,
            first: AtomRef::Index(0),
            second: AtomRef::Index(1),
            bond: BondDsl(BondAst::from_order(1)),
        })
    )]
    #[case::add_map_keyword(
        r##"{:bond {:add {:id :b1 :atoms [:c :nu] :type "2"}}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            keyword: Some("b1".into()),
            first: AtomRef::Keyword("c".into()),
            second: AtomRef::Keyword("nu".into()),
            bond: BondDsl(BondAst::from_order(2)),
        })
    )]
    #[case::remove_keyword("{:bond {:remove :b1}}", DeltaInput::BondRemove(BondRef::Keyword("b1".into())))]
    #[case::remove_index("{:bond {:remove 0}}", DeltaInput::BondRemove(BondRef::Index(0)))]
    #[case::modify(
        r##"{:bond {:modify [:b1 "2"]}}"##,
        DeltaInput::BondModify(
            BondRef::Keyword("b1".into()),
            BondUpdate { order: Some(ValueAst::Lit(2)), ..Default::default() },
        )
    )]
    #[case::modify_undetermined(
        r##"{:bond {:modify [:b1 "*#c*#u*#s*"]}}"##,
        DeltaInput::BondModify(
            BondRef::Keyword("b1".into()),
            BondUpdate {
                order: Some(ValueAst::Undetermined),
                charge: Some(ValueAst::Undetermined),
                spin: SpinStateUpdate {
                    unpaired: Some(ValueAst::Undetermined),
                    multiplicity: Some(ValueAst::Undetermined),
                },
                ..Default::default()
            },
        )
    )]
    fn test_read_delta_input_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(r##"{:dative-bond {:modify [:d1 "2"]}}"##, DeltaInput::DativeBondModify(DativeBondRef::Keyword("d1".into()), DativeBondUpdate { order: Some(ValueAst::Lit(2)), ..Default::default() }))]
    #[case::order_undetermined(r##"{:dative-bond {:modify [:d1 "*"]}}"##, DeltaInput::DativeBondModify(DativeBondRef::Keyword("d1".into()), DativeBondUpdate { order: Some(ValueAst::Undetermined), ..Default::default() }))]
    #[case::constraint_removal(r##"{:dative-bond {:modify [:d1 "#R(6)*"]}}"##, DeltaInput::DativeBondModify(DativeBondRef::Keyword("d1".into()), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() }))]
    fn test_parse_delta_input_dative_bond(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(r##"{:dative-bond {:modify [:d1 "2"]}}"##, DeltaInput::DativeBondModify(DativeBondRef::Keyword("d1".into()), DativeBondUpdate { order: Some(ValueAst::Lit(2)), ..Default::default() }))]
    #[case::order_undetermined(r##"{:dative-bond {:modify [:d1 "*"]}}"##, DeltaInput::DativeBondModify(DativeBondRef::Keyword("d1".into()), DativeBondUpdate { order: Some(ValueAst::Undetermined), ..Default::default() }))]
    #[case::constraint_removal(r##"{:dative-bond {:modify [:d1 "#R(6)*"]}}"##, DeltaInput::DativeBondModify(DativeBondRef::Keyword("d1".into()), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() }))]
    fn test_read_delta_input_dative_bond(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::spin_component(r##"{:aromatic-system {:modify [:a1 "#s1"]}}"##, DeltaInput::AromaticSystemModify(AromaticSystemRef::Keyword("a1".into()), AromaticSystemUpdate { spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined(r##"{:aromatic-system {:modify [:a1 "*#c*#e*"]}}"##, DeltaInput::AromaticSystemModify(AromaticSystemRef::Keyword("a1".into()), AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_parse_delta_input_aromatic_system(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::spin_component(r##"{:aromatic-system {:modify [:a1 "#s1"]}}"##, DeltaInput::AromaticSystemModify(AromaticSystemRef::Keyword("a1".into()), AromaticSystemUpdate { spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined(r##"{:aromatic-system {:modify [:a1 "*#c*#e*"]}}"##, DeltaInput::AromaticSystemModify(AromaticSystemRef::Keyword("a1".into()), AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_read_delta_input_aromatic_system(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::spin_component(r##"{:multicenter-bond {:modify [:m1 "#s1"]}}"##, DeltaInput::MulticenterBondModify(MulticenterBondRef::Keyword("m1".into()), MulticenterBondUpdate { spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined(r##"{:multicenter-bond {:modify [:m1 "*#c*#e*"]}}"##, DeltaInput::MulticenterBondModify(MulticenterBondRef::Keyword("m1".into()), MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_parse_delta_input_multicenter_bond(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::spin_component(r##"{:multicenter-bond {:modify [:m1 "#s1"]}}"##, DeltaInput::MulticenterBondModify(MulticenterBondRef::Keyword("m1".into()), MulticenterBondUpdate { spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined(r##"{:multicenter-bond {:modify [:m1 "*#c*#e*"]}}"##, DeltaInput::MulticenterBondModify(MulticenterBondRef::Keyword("m1".into()), MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_read_delta_input_multicenter_bond(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_undetermined(r##"{:noncovalent-bond {:modify [:n1 "*"]}}"##, DeltaInput::NoncovalentBondModify(NoncovalentBondRef::Keyword("n1".into()), NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Undetermined), ..Default::default() }))]
    #[case::constraint_removal(r##"{:noncovalent-bond {:modify [:n1 "#I*"]}}"##, DeltaInput::NoncovalentBondModify(NoncovalentBondRef::Keyword("n1".into()), NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(BooleanAst::Undetermined)), ..Default::default() }))]
    fn test_parse_delta_input_noncovalent_bond(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_undetermined(r##"{:noncovalent-bond {:modify [:n1 "*"]}}"##, DeltaInput::NoncovalentBondModify(NoncovalentBondRef::Keyword("n1".into()), NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Undetermined), ..Default::default() }))]
    #[case::constraint_removal(r##"{:noncovalent-bond {:modify [:n1 "#I*"]}}"##, DeltaInput::NoncovalentBondModify(NoncovalentBondRef::Keyword("n1".into()), NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(BooleanAst::Undetermined)), ..Default::default() }))]
    fn test_read_delta_input_noncovalent_bond(
        #[case] input: &str,
        #[case] expected: DeltaInput,
    ) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(r##"{:stereo-atom {:modify [:s1 ""]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate::default()))]
    #[case::undetermined(r##"{:stereo-atom {:modify [:s1 "*"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
    #[case::absolute(r##"{:stereo-atom {:modify [:s1 "Th1"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
    #[case::relative(r##"{:stereo-atom {:modify [:s1 "Th"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() }))]
    #[case::constraint_removal(r##"{:stereo-atom {:modify [:s1 "Th#g*"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)) }))]
    fn test_parse_delta_input_stereo_atom(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(parse_delta_input(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(r##"{:stereo-atom {:modify [:s1 ""]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate::default()))]
    #[case::undetermined(r##"{:stereo-atom {:modify [:s1 "*"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
    #[case::absolute(r##"{:stereo-atom {:modify [:s1 "Th1"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
    #[case::relative(r##"{:stereo-atom {:modify [:s1 "Th"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() }))]
    #[case::constraint_removal(r##"{:stereo-atom {:modify [:s1 "Th#g*"]}}"##, DeltaInput::StereoAtomModify(StereoAtomRef::Keyword("s1".into()), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)) }))]
    fn test_read_delta_input_stereo_atom(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(r##"{:stereo-bond {:modify [:t1 ""]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate::default()))]
    #[case::undetermined(r##"{:stereo-bond {:modify [:t1 "*"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
    #[case::absolute(r##"{:stereo-bond {:modify [:t1 "Ct1"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
    #[case::relative(r##"{:stereo-bond {:modify [:t1 "Ct"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() }))]
    #[case::constraint_removal(r##"{:stereo-bond {:modify [:t1 "Ct#g*"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)) }))]
    fn test_parse_delta_input_stereo_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(parse_delta_input(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(r##"{:stereo-bond {:modify [:t1 ""]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate::default()))]
    #[case::undetermined(r##"{:stereo-bond {:modify [:t1 "*"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
    #[case::absolute(r##"{:stereo-bond {:modify [:t1 "Ct1"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
    #[case::relative(r##"{:stereo-bond {:modify [:t1 "Ct"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() }))]
    #[case::constraint_removal(r##"{:stereo-bond {:modify [:t1 "Ct#g*"]}}"##, DeltaInput::StereoBondModify(StereoBondRef::Keyword("t1".into()), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)) }))]
    fn test_read_delta_input_stereo_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::add(
        "{:constraint {:add {:connected {}}}}",
        DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    #[case::remove(
        "{:constraint {:remove {:connected {}}}}",
        DeltaInput::ConstraintRemove(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    fn test_parse_delta_input_constraint(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            parse_delta_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::add(
        "{:constraint {:add {:connected {}}}}",
        DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    #[case::remove(
        "{:constraint {:remove {:connected {}}}}",
        DeltaInput::ConstraintRemove(ConstraintDsl::Molecule(MoleculeConstraintDsl::Connected {
            atoms: None,
        }))
    )]
    fn test_read_delta_input_constraint(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_parse_reaction_input() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:bond {:remove 0}} {:constraint {:add {:connected {}}}}]}"##;
        let expected = ReactionInput {
            lhs: MoleculeInput {
                atoms: vec![AtomEntryInput {
                    keyword: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))),
                }],
                ..Default::default()
            },
            atom_aliases: Vec::new(),
            deltas: vec![
                DeltaInput::AtomAdd(AtomEntryInput {
                    keyword: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))),
                }),
                DeltaInput::BondRemove(BondRef::Index(0)),
                DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(
                    MoleculeConstraintDsl::Connected { atoms: None },
                )),
            ],
        };
        assert_eq!(
            parse_reaction_input(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_read_reaction_input() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:bond {:remove 0}} {:constraint {:add {:connected {}}}}]}"##;
        let expected = ReactionInput {
            lhs: MoleculeInput {
                atoms: vec![AtomEntryInput {
                    keyword: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::C)))),
                }],
                ..Default::default()
            },
            atom_aliases: Vec::new(),
            deltas: vec![
                DeltaInput::AtomAdd(AtomEntryInput {
                    keyword: None,
                    spec: AtomSpecInput::Bare(Box::new(AtomDsl(AtomAst::from_element(Element::O)))),
                }),
                DeltaInput::BondRemove(BondRef::Index(0)),
                DeltaInput::ConstraintAdd(ConstraintDsl::Molecule(
                    MoleculeConstraintDsl::Connected { atoms: None },
                )),
            ],
        };
        assert_eq!(
            read_reaction_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast() {
        let input = r##"{:lhs {:atoms ["C"]} :atom-aliases [:me "C#h3"] :deltas [{:atom {:add [:nu :me]}} {:atom {:add "O"}}]}"##;
        let (ast, meta) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: {
                        let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                        a.implicit_hydrogens = ValueAst::Lit(3);
                        a
                    },
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::O),
                }),
            ]),
        );
        assert_eq!(meta.atom_keyword(AtomId(1)), Some("nu"));
        assert_eq!(meta.atom_keyword(AtomId(2)), None);
        assert!(meta.has_atom_alias("me"));
    }

    #[rstest]
    fn test_reaction_input_into_ast_alias_union() {
        // `:lo` is an lhs alias, `:hi` a reaction alias; both `:add` resolve (union),
        // but each set stays in its own metadata field for independent round-trip.
        let input = r##"{:lhs {:atoms ["C"] :atom-aliases [:lo "N"]} :atom-aliases [:hi "C#h3"] :deltas [{:atom {:add :lo}} {:atom {:add :hi}}]}"##;
        let (ast, meta) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: {
                        let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                        a.implicit_hydrogens = ValueAst::Lit(3);
                        a
                    },
                }),
            ]),
        );
        assert_eq!(meta.lhs().atom_aliases_len(), 1);
        assert!(meta.lhs().has_atom_alias("lo"));
        assert_eq!(meta.atom_aliases_len(), 1);
        assert!(meta.has_atom_alias("hi"));
    }

    #[rstest]
    fn test_reaction_input_into_ast_atom_remove() {
        let input = r##"{:lhs {:atoms [[:br "Br"] "C"]} :deltas [{:atom {:remove :br}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::Br),
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_atom_remove_error() {
        // Adding then removing through the same keyword is prohibited:
        // recover-from-lhs cannot reach an added atom.
        let input =
            r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:x "O"]}} {:atom {:remove :x}}]}"##;
        let err = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::DeltaTargetAdded {
                action: "remove",
                kind: "atom",
                index: 1
            }
        );
    }

    // A delta keyword must be disjoint from every keyword already bound, whether by an lhs entity
    // or an earlier delta.
    #[rstest]
    #[case::collides_with_lhs(
        r##"{:lhs {:atoms [[:a "C"]]} :deltas [{:atom {:add [:a "O"]}}]}"##,
        ParseError::DuplicateKeyword("a".to_string()),
    )]
    #[case::collides_with_prior_delta(
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:a "O"]}} {:atom {:add [:a "N"]}}]}"##,
        ParseError::DuplicateKeyword("a".to_string()),
    )]
    fn test_reaction_input_into_ast_duplicate_keyword_error(
        #[case] input: &str,
        #[case] expected: ParseError,
    ) {
        assert_eq!(
            parse_reaction_input(&read_string(input).unwrap())
                .unwrap()
                .into_ast()
                .unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::determined(
        r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##,
        ValueAst::Lit(-1),
    )]
    #[case::undetermined(
        r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c*"]}}]}"##,
        ValueAst::Undetermined
    )]
    fn test_reaction_input_into_ast_atom_modify(#[case] input: &str, #[case] new: ValueAst) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: ValueAst::Lit(0),
                    new,
                },
            })]),
        );
    }

    #[rstest]
    #[case::multiplicity(
        r##"{:lhs {:atoms [[:c "C#u2#s3"]]} :deltas [{:atom {:modify [:c "#s1"]}}]}"##,
        SpinStateAst::from((2_u8, 1_u8)),
    )]
    fn test_reaction_input_into_ast_atom_modify_spin(
        #[case] input: &str,
        #[case] new: SpinStateAst,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Spin {
                    old: SpinStateAst::from((2_u8, 3_u8)),
                    new,
                },
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_atom_modify_remove_constraint() {
        // An undetermined constraint in an AtomUpdate removes that keyed constraint.
        let input = r##"{:lhs {:atoms [[:me "C#v4"]]} :deltas [{:atom {:modify [:me "#v*"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: Some(AtomConstraintAst::valence(4)),
                new: None,
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_add() {
        // Bond-add attaches to a same-reaction atom: :o is AtomId(1), bond endpoints [0, :o].
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:bond {:add [0 :o "1"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_remove() {
        let input = r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:remove :b1}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondAst::from_order(1),
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_remove_error() {
        // Adding then removing the same bond is prohibited — recover-from-lhs cannot reach it.
        let input =
            r##"{:lhs {:atoms ["C" "O"]} :deltas [{:bond {:add [0 1 "1"]}} {:bond {:remove 0}}]}"##;
        let err = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::DeltaTargetAdded {
                action: "remove",
                kind: "bond",
                index: 0
            }
        );
    }

    #[rstest]
    #[case::determined(
        r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "2"]}}]}"##,
        ValueAst::Lit(2),
    )]
    #[case::undetermined(
        r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "*"]}}]}"##,
        ValueAst::Undetermined,
    )]
    fn test_reaction_input_into_ast_bond_modify(#[case] input: &str, #[case] new: ValueAst) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new,
                },
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_modify_spin() {
        let input = r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1#u2#s3"}]} :deltas [{:bond {:modify [:b1 "#s1"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Spin {
                    old: SpinStateAst::from((2_u8, 3_u8)),
                    new: SpinStateAst::from((2_u8, 1_u8)),
                },
            })]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_bond_modify_constraint() {
        let input = r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1#R(6)"}]} :deltas [{:bond {:modify [:b1 "#R(6)*"]}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyConstraint {
                id: BondId(0),
                old: Some(BondConstraintAst::ring_membership(
                    RingScope::Size(6),
                    ValueAst::Lit(1),
                )),
                new: None,
            })]),
        );
    }

    #[rstest]
    #[case::determined(
        r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1"}]} :deltas [{:dative-bond {:modify [:d1 "2"]}}]}"##,
        ValueAst::Lit(2),
    )]
    #[case::undetermined(
        r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1"}]} :deltas [{:dative-bond {:modify [:d1 "*"]}}]}"##,
        ValueAst::Undetermined,
    )]
    fn test_reaction_input_into_ast_dative_bond_modify(#[case] input: &str, #[case] new: ValueAst) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::DativeBond(DativeBondDelta::ModifyField {
                id: DativeBondId(0),
                change: DativeBondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new,
                },
            })]),
        );
    }

    #[rstest]
    #[case::ring_size(
        r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:id :d1 :donors [0] :acceptor 1 :type "1#R(6)"}]} :deltas [{:dative-bond {:modify [:d1 "#R(6)*"]}}]}"##,
        DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Lit(1)),
    )]
    fn test_reaction_input_into_ast_dative_bond_modify_constraint(
        #[case] input: &str,
        #[case] old: DativeBondConstraintAst,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::DativeBond(DativeBondDelta::ModifyConstraint {
                id: DativeBondId(0),
                old: Some(old),
                new: None,
            },)]),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons_undetermined(
        r##"{:lhs {:atoms ["C" "C" "C"] :aromatic-systems [{:id :a1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:aromatic-system {:modify [:a1 "*"]}}]}"##,
        vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyField { id: AromaticSystemId(0), change: AromaticSystemFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1, 1]), new: ElectronCountsAst::Undetermined } })],
    )]
    #[case::charge_undetermined(
        r##"{:lhs {:atoms ["C" "C" "C"] :aromatic-systems [{:id :a1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:aromatic-system {:modify [:a1 "#c*"]}}]}"##,
        vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyField { id: AromaticSystemId(0), change: AromaticSystemFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined } })],
    )]
    #[case::spin_component(
        r##"{:lhs {:atoms ["C" "C" "C"] :aromatic-systems [{:id :a1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:aromatic-system {:modify [:a1 "#s1"]}}]}"##,
        vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyField { id: AromaticSystemId(0), change: AromaticSystemFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) } })],
    )]
    #[case::constraint_removal(
        r##"{:lhs {:atoms ["C" "C" "C"] :aromatic-systems [{:id :a1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:aromatic-system {:modify [:a1 "#e*"]}}]}"##,
        vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyConstraint { id: AromaticSystemId(0), old: Some(AromaticSystemConstraintAst::electron_count(6_i64)), new: None })],
    )]
    fn test_reaction_input_into_ast_aromatic_system_modify(
        #[case] input: &str,
        #[case] expected: Vec<Delta>,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(ast.deltas, Deltas::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons_undetermined(
        r##"{:lhs {:atoms ["C" "C" "C"] :multicenter-bonds [{:id :m1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:multicenter-bond {:modify [:m1 "*"]}}]}"##,
        vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyField { id: MulticenterBondId(0), change: MulticenterBondFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1, 1]), new: ElectronCountsAst::Undetermined } })],
    )]
    #[case::charge_undetermined(
        r##"{:lhs {:atoms ["C" "C" "C"] :multicenter-bonds [{:id :m1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:multicenter-bond {:modify [:m1 "#c*"]}}]}"##,
        vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyField { id: MulticenterBondId(0), change: MulticenterBondFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined } })],
    )]
    #[case::spin_component(
        r##"{:lhs {:atoms ["C" "C" "C"] :multicenter-bonds [{:id :m1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:multicenter-bond {:modify [:m1 "#s1"]}}]}"##,
        vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyField { id: MulticenterBondId(0), change: MulticenterBondFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) } })],
    )]
    #[case::constraint_removal(
        r##"{:lhs {:atoms ["C" "C" "C"] :multicenter-bonds [{:id :m1 :atoms [0 1 2] :type "[1,1,1]#c0#u2#s3#e6"}]} :deltas [{:multicenter-bond {:modify [:m1 "#e*"]}}]}"##,
        vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyConstraint { id: MulticenterBondId(0), old: Some(MulticenterBondConstraintAst::electron_count(6_i64)), new: None })],
    )]
    fn test_reaction_input_into_ast_multicenter_bond_modify(
        #[case] input: &str,
        #[case] expected: Vec<Delta>,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(ast.deltas, Deltas::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_undetermined(
        r##"{:lhs {:atoms ["N" "H"] :noncovalent-bonds [{:id :n1 :atoms [0 1] :type "Hbd"}]} :deltas [{:noncovalent-bond {:modify [:n1 "*"]}}]}"##,
        vec![Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField { id: NoncovalentBondId(0), change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindAst::Undetermined } })],
    )]
    #[case::constraint_removal(
        r##"{:lhs {:atoms ["N" "H"] :noncovalent-bonds [{:id :n1 :atoms [0 1] :type "Hbd#I"}]} :deltas [{:noncovalent-bond {:modify [:n1 "#I*"]}}]}"##,
        vec![Delta::NoncovalentBond(NoncovalentBondDelta::ModifyConstraint { id: NoncovalentBondId(0), old: Some(NoncovalentBondConstraintAst::intramolecular(true)), new: None })],
    )]
    fn test_reaction_input_into_ast_noncovalent_bond_modify(
        #[case] input: &str,
        #[case] expected: Vec<Delta>,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(ast.deltas, Deltas::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 3 4] :type "Th0"}]} :deltas [{:stereo-atom {:modify [:s1 "Th1"]}}]}"##,
        vec![Delta::StereoAtom(StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 1_u32) } })],
    )]
    #[case::undetermined(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 3 4] :type "Th0"}]} :deltas [{:stereo-atom {:modify [:s1 "*"]}}]}"##,
        vec![Delta::StereoAtom(StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationAst::Undetermined } })],
    )]
    #[case::constraint_removal(
        r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :stereo-atoms [{:id :s1 :site 0 :ligands [1 2 3 4] :type "Th0#g/"}]} :deltas [{:stereo-atom {:modify [:s1 "Th#g*"]}}]}"##,
        vec![Delta::StereoAtom(StereoAtomDelta::ModifyConstraint { id: StereoAtomId(0), kind: Some(StereoKind::Tetrahedral), old: Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))), new: None })],
    )]
    fn test_reaction_input_into_ast_stereo_atom_modify(
        #[case] input: &str,
        #[case] expected: Vec<Delta>,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(ast.deltas, Deltas::from_iter(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(
        r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :t1 :site 1 :ligands [0 3] :type "Ct0"}]} :deltas [{:stereo-bond {:modify [:t1 "Ct1"]}}]}"##,
        vec![Delta::StereoBond(StereoBondDelta::ModifyField { id: StereoBondId(0), change: StereoBondFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationAst::kinded(StereoKind::CisTrans, 1_u32) } })],
    )]
    #[case::undetermined(
        r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :t1 :site 1 :ligands [0 3] :type "Ct0"}]} :deltas [{:stereo-bond {:modify [:t1 "*"]}}]}"##,
        vec![Delta::StereoBond(StereoBondDelta::ModifyField { id: StereoBondId(0), change: StereoBondFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationAst::Undetermined } })],
    )]
    #[case::constraint_removal(
        r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:id :t1 :site 1 :ligands [0 3] :type "Ct0#g/"}]} :deltas [{:stereo-bond {:modify [:t1 "Ct#g*"]}}]}"##,
        vec![Delta::StereoBond(StereoBondDelta::ModifyConstraint { id: StereoBondId(0), kind: Some(StereoKind::CisTrans), old: Some(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))), new: None })],
    )]
    fn test_reaction_input_into_ast_stereo_bond_modify(
        #[case] input: &str,
        #[case] expected: Vec<Delta>,
    ) {
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(ast.deltas, Deltas::from_iter(expected));
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_add() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:add {:connected {}}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None },)
            ))]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_remove() {
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:constraint {:remove {:connected {}}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Remove(
                Constraint::Molecule(MoleculeConstraint::Connected { atoms: None },)
            ))]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_added_atom_ref() {
        // The constraint ref :o names an atom added in the same reaction (AtomId(1)), resolved
        // against the unified namespace.
        let input = r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:constraint {:add {:atom [:o {:valence 2}]}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Constraint(ConstraintDelta::Add(Constraint::Atom(
                    AtomId(1),
                    AtomConstraintAst::Valence(ValueAst::Lit(2)),
                ))),
            ]),
        );
    }

    #[rstest]
    fn test_reaction_input_into_ast_constraint_structural_bond_ref() {
        // A structural bond ref ({:atoms [0 1]}) names the lhs bond by its endpoints, resolved
        // against the namespace's participant lookup.
        let input = r##"{:lhs {:atoms ["C" "C"] :bonds [[0 1 "1"]]} :deltas [{:constraint {:add {:bond [{:atoms [0 1]} {:aromatic true}]}}}]}"##;
        let (ast, _) = parse_reaction_input(&read_string(input).unwrap())
            .unwrap()
            .into_ast()
            .unwrap();
        assert_eq!(
            ast.deltas,
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(Constraint::Bond(
                BondId(0),
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            )))]),
        );
    }

    #[rstest]
    #[case::atom_modify(
        r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##,
        ReactionAst {
            lhs: MoleculeAst::from_edn_str(r##"{:atoms [[:br "Br#c0"]]}"##).unwrap(),
            deltas: Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) },
            })]),
        }
    )]
    #[case::atom_add_bond_add(
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:bond {:add [0 :o "1"]}}]}"##,
        ReactionAst {
            lhs: MoleculeAst::from_edn_str(r##"{:atoms ["C"]}"##).unwrap(),
            deltas: Deltas::from_iter([
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        }
    )]
    fn test_reaction_dsl_from_edn(#[case] input: &str, #[case] expected: ReactionAst) {
        let dsl = ReactionDsl::from_edn(&read_string(input).unwrap()).unwrap();
        assert_eq!(dsl.ast(), &expected);
    }

    #[rstest]
    #[case::atom_modify(
        r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##
    )]
    #[case::bond_modify_and_constraint(
        r##"{:lhs {:atoms ["C" "O"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "2"]}} {:constraint {:add {:connected {}}}}]}"##
    )]
    #[case::atom_add_bond_add(
        r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:bond {:add [0 :o "1"]}}]}"##
    )]
    fn test_reaction_dsl_from_edn_str_from_edn_parity(#[case] input: &str) {
        let via_tree = ReactionDsl::from_edn(&read_string(input).unwrap()).unwrap();
        let via_stream = ReactionDsl::from_edn_str(input).unwrap();
        assert_eq!(via_tree, via_stream);
    }

    #[rstest]
    fn test_reaction_dsl_from_str() {
        let input = r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##;
        let dsl: ReactionDsl = input.parse().unwrap();
        assert_eq!(dsl, ReactionDsl::from_edn_str(input).unwrap());
    }

    /// Shared render metadata for existing and created test entities.
    #[fixture]
    fn meta() -> ReactionMetadata {
        let mut lhs = MoleculeMetadata::new();
        lhs.set_keyword(Entity::Atom(AtomId(0)), "br").unwrap();
        lhs.set_keyword(Entity::Atom(AtomId(1)), "c").unwrap();
        lhs.set_keyword(Entity::Bond(BondId(0)), "b1").unwrap();
        lhs.set_keyword(Entity::Bond(BondId(1)), "bx").unwrap();
        lhs.set_keyword(Entity::DativeBond(DativeBondId(0)), "d1")
            .unwrap();
        lhs.set_keyword(Entity::AromaticSystem(AromaticSystemId(0)), "a1")
            .unwrap();
        lhs.set_keyword(Entity::MulticenterBond(MulticenterBondId(0)), "m1")
            .unwrap();
        lhs.set_keyword(Entity::NoncovalentBond(NoncovalentBondId(0)), "n1")
            .unwrap();
        lhs.set_keyword(Entity::StereoAtom(StereoAtomId(0)), "s1")
            .unwrap();
        lhs.set_keyword(Entity::StereoBond(StereoBondId(0)), "t1")
            .unwrap();

        let mut metadata = ReactionMetadata::from(lhs);
        metadata.set_atom_keyword(AtomId(2), "n");
        metadata.set_bond_keyword(BondId(2), "b2");
        metadata
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::add(vec![Delta::Atom(AtomDelta::Add { id: AtomId(2), ast: AtomAst::from_element(Element::N) })], r##"{:atom {:add [:n "N"]}}"##)]
    #[case::remove(vec![Delta::Atom(AtomDelta::Remove { id: AtomId(1), ast: AtomAst::from_element(Element::C) })], "{:atom {:remove :c}}")]
    #[case::modify_field(vec![Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) } })], r##"{:atom {:modify [:br "#c-"]}}"##)]
    #[case::modify_field_undetermined(vec![Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined } })], r##"{:atom {:modify [:br "#c*"]}}"##)]
    #[case::modify_spin(vec![Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) } })], r##"{:atom {:modify [:br "#u2#s"]}}"##)]
    #[case::modify_set_constraint(vec![Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: None, new: Some(AtomConstraintAst::valence(4_i64)) })], r##"{:atom {:modify [:br "#v4"]}}"##)]
    #[case::modify_remove_constraint(vec![Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: Some(AtomConstraintAst::valence(4_i64)), new: None })], r##"{:atom {:modify [:br "#v*"]}}"##)]
    #[case::modify_coalesced(vec![Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(-1) } }), Delta::Atom(AtomDelta::ModifyConstraint { id: AtomId(0), old: Some(AtomConstraintAst::valence(4_i64)), new: None })], r##"{:atom {:modify [:br "#c-#v*"]}}"##)]
    fn test_render_deltas_atom(meta: ReactionMetadata, #[case] deltas: Vec<Delta>, #[case] expected: &str) {
        assert_eq!(render_deltas(&Deltas::from_iter(deltas), &meta), vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::add(vec![Delta::Bond(BondDelta::Add { id: BondId(2), atoms: [AtomId(1), AtomId(2)], ast: BondAst::from_order(1) })], "{:bond {:add {:id :b2 :atoms [:c :n] :type :single}}}")]
    #[case::remove(vec![Delta::Bond(BondDelta::Remove { id: BondId(1), atoms: [AtomId(0), AtomId(1)], ast: BondAst::from_order(1) })], "{:bond {:remove :bx}}")]
    #[case::modify_field(vec![Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) } })], r##"{:bond {:modify [:b1 "2"]}}"##)]
    #[case::modify_field_undetermined(vec![Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Undetermined } })], r##"{:bond {:modify [:b1 "*"]}}"##)]
    #[case::modify_charge_undetermined(vec![Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined } })], r##"{:bond {:modify [:b1 "#c*"]}}"##)]
    #[case::modify_spin(vec![Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) } })], r##"{:bond {:modify [:b1 "#u2#s"]}}"##)]
    #[case::modify_constraint(vec![Delta::Bond(BondDelta::ModifyConstraint { id: BondId(0), old: None, new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))) })], r##"{:bond {:modify [:b1 "#a"]}}"##)]
    #[case::modify_remove_constraint(vec![Delta::Bond(BondDelta::ModifyConstraint { id: BondId(0), old: Some(BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Lit(1))), new: None })], r##"{:bond {:modify [:b1 "#R(6)*"]}}"##)]
    fn test_render_deltas_bond(meta: ReactionMetadata, #[case] deltas: Vec<Delta>, #[case] expected: &str) {
        assert_eq!(render_deltas(&Deltas::from_iter(deltas), &meta), vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(vec![Delta::DativeBond(DativeBondDelta::ModifyField { id: DativeBondId(0), change: DativeBondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) } })], r##"{:dative-bond {:modify [:d1 "2"]}}"##)]
    #[case::order_undetermined(vec![Delta::DativeBond(DativeBondDelta::ModifyField { id: DativeBondId(0), change: DativeBondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Undetermined } })], r##"{:dative-bond {:modify [:d1 "*"]}}"##)]
    #[case::constraint(vec![Delta::DativeBond(DativeBondDelta::ModifyConstraint { id: DativeBondId(0), old: None, new: Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) })], r##"{:dative-bond {:modify [:d1 "#a"]}}"##)]
    #[case::constraint_removal(vec![Delta::DativeBond(DativeBondDelta::ModifyConstraint { id: DativeBondId(0), old: Some(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Lit(1))), new: None })], r##"{:dative-bond {:modify [:d1 "#R(6)*"]}}"##)]
    fn test_render_deltas_dative_bond(
        meta: ReactionMetadata,
        #[case] deltas: Vec<Delta>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_deltas(&Deltas::from_iter(deltas), &meta),
            vec![read_string(expected).unwrap()],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons_undetermined(vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyField { id: AromaticSystemId(0), change: AromaticSystemFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1, 1]), new: ElectronCountsAst::Undetermined } })], r##"{:aromatic-system {:modify [:a1 "*"]}}"##)]
    #[case::charge_undetermined(vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyField { id: AromaticSystemId(0), change: AromaticSystemFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined } })], r##"{:aromatic-system {:modify [:a1 "#c*"]}}"##)]
    #[case::spin(vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyField { id: AromaticSystemId(0), change: AromaticSystemFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) } })], r##"{:aromatic-system {:modify [:a1 "#u2#s"]}}"##)]
    #[case::constraint_removal(vec![Delta::AromaticSystem(AromaticSystemDelta::ModifyConstraint { id: AromaticSystemId(0), old: Some(AromaticSystemConstraintAst::electron_count(6_i64)), new: None })], r##"{:aromatic-system {:modify [:a1 "#e*"]}}"##)]
    fn test_render_deltas_aromatic_system(
        meta: ReactionMetadata,
        #[case] deltas: Vec<Delta>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_deltas(&Deltas::from_iter(deltas), &meta),
            vec![read_string(expected).unwrap()],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons_undetermined(vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyField { id: MulticenterBondId(0), change: MulticenterBondFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1, 1]), new: ElectronCountsAst::Undetermined } })], r##"{:multicenter-bond {:modify [:m1 "*"]}}"##)]
    #[case::charge_undetermined(vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyField { id: MulticenterBondId(0), change: MulticenterBondFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined } })], r##"{:multicenter-bond {:modify [:m1 "#c*"]}}"##)]
    #[case::spin(vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyField { id: MulticenterBondId(0), change: MulticenterBondFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) } })], r##"{:multicenter-bond {:modify [:m1 "#u2#s"]}}"##)]
    #[case::constraint_removal(vec![Delta::MulticenterBond(MulticenterBondDelta::ModifyConstraint { id: MulticenterBondId(0), old: Some(MulticenterBondConstraintAst::electron_count(6_i64)), new: None })], r##"{:multicenter-bond {:modify [:m1 "#e*"]}}"##)]
    fn test_render_deltas_multicenter_bond(
        meta: ReactionMetadata,
        #[case] deltas: Vec<Delta>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_deltas(&Deltas::from_iter(deltas), &meta),
            vec![read_string(expected).unwrap()],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_undetermined(vec![Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField { id: NoncovalentBondId(0), change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindAst::Undetermined } })], r##"{:noncovalent-bond {:modify [:n1 "*"]}}"##)]
    #[case::kind(vec![Delta::NoncovalentBond(NoncovalentBondDelta::ModifyField { id: NoncovalentBondId(0), change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic) } })], r##"{:noncovalent-bond {:modify [:n1 "Ion"]}}"##)]
    #[case::constraint_removal(vec![Delta::NoncovalentBond(NoncovalentBondDelta::ModifyConstraint { id: NoncovalentBondId(0), old: Some(NoncovalentBondConstraintAst::intramolecular(true)), new: None })], r##"{:noncovalent-bond {:modify [:n1 "#I*"]}}"##)]
    fn test_render_deltas_noncovalent_bond(
        meta: ReactionMetadata,
        #[case] deltas: Vec<Delta>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_deltas(&Deltas::from_iter(deltas), &meta),
            vec![read_string(expected).unwrap()],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(vec![Delta::StereoAtom(StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 1_u32) } })], r##"{:stereo-atom {:modify [:s1 "Th1"]}}"##)]
    #[case::undetermined(vec![Delta::StereoAtom(StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationAst::Undetermined } })], r##"{:stereo-atom {:modify [:s1 "*"]}}"##)]
    #[case::constraint_removal(vec![Delta::StereoAtom(StereoAtomDelta::ModifyConstraint { id: StereoAtomId(0), kind: Some(StereoKind::Tetrahedral), old: Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))), new: None })], r##"{:stereo-atom {:modify [:s1 "Th#g*"]}}"##)]
    fn test_render_deltas_stereo_atom(
        meta: ReactionMetadata,
        #[case] deltas: Vec<Delta>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_deltas(&Deltas::from_iter(deltas), &meta),
            vec![read_string(expected).unwrap()],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(vec![Delta::StereoBond(StereoBondDelta::ModifyField { id: StereoBondId(0), change: StereoBondFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationAst::kinded(StereoKind::CisTrans, 1_u32) } })], r##"{:stereo-bond {:modify [:t1 "Ct1"]}}"##)]
    #[case::undetermined(vec![Delta::StereoBond(StereoBondDelta::ModifyField { id: StereoBondId(0), change: StereoBondFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationAst::Undetermined } })], r##"{:stereo-bond {:modify [:t1 "*"]}}"##)]
    #[case::constraint_removal(vec![Delta::StereoBond(StereoBondDelta::ModifyConstraint { id: StereoBondId(0), kind: Some(StereoKind::CisTrans), old: Some(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))), new: None })], r##"{:stereo-bond {:modify [:t1 "Ct#g*"]}}"##)]
    fn test_render_deltas_stereo_bond(
        meta: ReactionMetadata,
        #[case] deltas: Vec<Delta>,
        #[case] expected: &str,
    ) {
        assert_eq!(
            render_deltas(&Deltas::from_iter(deltas), &meta),
            vec![read_string(expected).unwrap()],
        );
    }

    #[rstest]
    fn test_render_deltas_stereo_bond_configuration_clear_with_constraint_removal(
        meta: ReactionMetadata,
    ) {
        let deltas = Deltas::from_iter([
            Delta::StereoBond(StereoBondDelta::ModifyField {
                id: StereoBondId(0),
                change: StereoBondFieldChange::Configuration {
                    old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32),
                    new: StereoConfigurationAst::Undetermined,
                },
            }),
            Delta::StereoBond(StereoBondDelta::ModifyConstraint {
                id: StereoBondId(0),
                kind: Some(StereoKind::CisTrans),
                old: Some(StereoBondConstraintAst::Stereogenicity(
                    StereogenicityAst::Lit(Stereogenicity::Stereogenic),
                )),
                new: None,
            }),
        ]);
        assert_eq!(
            render_deltas(&deltas, &meta),
            vec![
                read_string(r##"{:stereo-bond {:modify [:t1 "*"]}}"##).unwrap(),
                read_string(r##"{:stereo-bond {:modify [:t1 "Ct#g*"]}}"##).unwrap(),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::add_molecule(vec![Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })))], "{:constraint {:add {:connected {}}}}")]
    #[case::add_entity_leaf(vec![Delta::Constraint(ConstraintDelta::Add(Constraint::Atom(AtomId(2), AtomConstraintAst::Valence(ValueAst::Lit(2)))))], "{:constraint {:add {:atom [:n {:valence 2}]}}}")]
    #[case::remove(vec![Delta::Constraint(ConstraintDelta::Remove(Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })))], "{:constraint {:remove {:connected {}}}}")]
    fn test_render_deltas_constraint(meta: ReactionMetadata, #[case] deltas: Vec<Delta>, #[case] expected: &str) {
        assert_eq!(render_deltas(&Deltas::from_iter(deltas), &meta), vec![read_string(expected).unwrap()]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c-1"]}}]}"##)]
    #[case::modify_undetermined(r##"{:lhs {:atoms [[:br "Br#c0"]]} :deltas [{:atom {:modify [:br "#c*"]}}]}"##)]
    #[case::modify_spin_component(r##"{:lhs {:atoms [[:c "C#u2#s3"]]} :deltas [{:atom {:modify [:c "#s1"]}}]}"##)]
    #[case::bond_modify_undetermined(r##"{:lhs {:atoms ["C" "N"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "*"]}}]}"##)]
    #[case::bond_modify_spin_component(r##"{:lhs {:atoms ["C" "N"] :bonds [{:id :b1 :atoms [0 1] :type "1#u2#s3"}]} :deltas [{:bond {:modify [:b1 "#s1"]}}]}"##)]
    #[case::bond_modify_constraint_removal(r##"{:lhs {:atoms ["C" "N"] :bonds [{:id :b1 :atoms [0 1] :type "1#R(6)"}]} :deltas [{:bond {:modify [:b1 "#R(6)*"]}}]}"##)]
    #[case::reaction_alias(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add :nu}}] :atom-aliases [:nu "O#h1#c-1"] }"##)]
    #[case::molecule_constraint(r##"{:lhs {:atoms ["C" "N"] :bonds [{:id :b1 :atoms [0 1] :type "1"}]} :deltas [{:bond {:modify [:b1 "2"]}} {:constraint {:add {:connected {}}}}]}"##)]
    #[case::entity_leaf_constraint(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add [:o "O"]}} {:constraint {:add {:atom [:o {:valence 2}]}}}]}"##)]
    #[case::dative_add(r##"{:lhs {:atoms ["C" "N"]} :deltas [{:dative-bond {:add {:donors [0] :acceptor 1 :type "1#R"}}}]}"##)]
    #[case::dative_remove(r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}]} :deltas [{:dative-bond {:remove 0}}]}"##)]
    #[case::dative_modify(r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R"}]} :deltas [{:dative-bond {:modify [0 "2"]}}]}"##)]
    #[case::dative_modify_undetermined(r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]} :deltas [{:dative-bond {:modify [0 "*"]}}]}"##)]
    #[case::dative_ring_constraint_removal(r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1#R(6)"}]} :deltas [{:dative-bond {:modify [0 "#R(6)*"]}}]}"##)]
    #[case::aromatic_add(r##"{:lhs {:atoms ["C" "C" "C" "C" "C" "C"]} :deltas [{:aromatic-system {:add {:atoms [0 1 2 3 4 5] :type "*#e6"}}}]}"##)]
    #[case::aromatic_remove(r##"{:lhs {:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "*#e6"}]} :deltas [{:aromatic-system {:remove 0}}]}"##)]
    #[case::multicenter_add(r##"{:lhs {:atoms ["C" "C"]} :deltas [{:multicenter-bond {:add {:atoms [0 1] :type "*#e2"}}}]}"##)]
    #[case::multicenter_remove(r##"{:lhs {:atoms ["C" "C"] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]} :deltas [{:multicenter-bond {:remove 0}}]}"##)]
    #[case::noncovalent_add(r##"{:lhs {:atoms ["N" "H"]} :deltas [{:noncovalent-bond {:add {:atoms [0 1] :type "Hbd"}}}]}"##)]
    #[case::noncovalent_remove(r##"{:lhs {:atoms ["N" "H"] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]} :deltas [{:noncovalent-bond {:remove 0}}]}"##)]
    #[case::noncovalent_modify_undetermined(r##"{:lhs {:atoms ["N" "H"] :noncovalent-bonds [{:atoms [0 1] :type "Hbd"}]} :deltas [{:noncovalent-bond {:modify [0 "*"]}}]}"##)]
    #[case::noncovalent_constraint_removal(r##"{:lhs {:atoms ["N" "H"] :noncovalent-bonds [{:atoms [0 1] :type "Hbd#I"}]} :deltas [{:noncovalent-bond {:modify [0 "#I*"]}}]}"##)]
    #[case::stereo_atom_add(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]} :deltas [{:stereo-atom {:add {:site 0 :ligands [1 2 3 4] :type "Th1"}}}]}"##)]
    #[case::stereo_atom_remove(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:remove 0}}]}"##)]
    #[case::stereo_atom_modify(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:modify [0 "Th2"]}}]}"##)]
    #[case::stereo_atom_modify_undetermined(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:modify [0 "*"]}}]}"##)]
    #[case::stereo_atom_swap(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:swap [0 :tetrahedral]}}]}"##)]
    #[case::stereo_atom_mirror(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:mirror [0 :tetrahedral]}}]}"##)]
    #[case::stereo_atom_apply(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]} :deltas [{:stereo-atom {:apply [0 :tetrahedral "(0,1)"]}}]}"##)]
    #[case::stereo_bond_add(r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]} :deltas [{:stereo-bond {:add {:site 1 :ligands [0 3] :type "Ct1"}}}]}"##)]
    #[case::stereo_bond_remove(r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]} :deltas [{:stereo-bond {:remove 0}}]}"##)]
    #[case::stereo_bond_modify(r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct0"}]} :deltas [{:stereo-bond {:modify [0 "Ct1"]}}]}"##)]
    #[case::stereo_bond_modify_undetermined(r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct0"}]} :deltas [{:stereo-bond {:modify [0 "*"]}}]}"##)]
    #[case::stereo_bond_stereogenicity_removal(r##"{:lhs {:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct0#g/"}]} :deltas [{:stereo-bond {:modify [0 "Ct#g*"]}}]}"##)]
    #[case::dative_constraint_removal(r##"{:lhs {:atoms ["C" "N"] :dative-bonds [{:donors [0] :acceptor 1 :type "1#a"}]} :deltas [{:dative-bond {:modify [0 "#a*"]}}]}"##)]
    #[case::aromatic_modify_undetermined(r##"{:lhs {:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]#c0#u2#s3#e6"}]} :deltas [{:aromatic-system {:modify [0 "*"]}}]}"##)]
    #[case::aromatic_modify_spin_component(r##"{:lhs {:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]#c0#u2#s3#e6"}]} :deltas [{:aromatic-system {:modify [0 "#s1"]}}]}"##)]
    #[case::aromatic_constraint_removal(r##"{:lhs {:atoms ["C" "C" "C" "C" "C" "C"] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "*#e6"}]} :deltas [{:aromatic-system {:modify [0 "#e*"]}}]}"##)]
    #[case::multicenter_modify_undetermined(r##"{:lhs {:atoms ["C" "C"] :multicenter-bonds [{:atoms [0 1] :type "[1,1]#c0#u2#s3#e2"}]} :deltas [{:multicenter-bond {:modify [0 "*"]}}]}"##)]
    #[case::multicenter_modify_spin_component(r##"{:lhs {:atoms ["C" "C"] :multicenter-bonds [{:atoms [0 1] :type "[1,1]#c0#u2#s3#e2"}]} :deltas [{:multicenter-bond {:modify [0 "#s1"]}}]}"##)]
    #[case::multicenter_constraint_removal(r##"{:lhs {:atoms ["C" "C"] :multicenter-bonds [{:atoms [0 1] :type "*#e2"}]} :deltas [{:multicenter-bond {:modify [0 "#e*"]}}]}"##)]
    #[case::stereo_atom_topicity_removal(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1#o(0,1)="}]} :deltas [{:stereo-atom {:modify [0 "Th#o(0,1)*"]}}]}"##)]
    #[case::stereo_atom_topicity_change(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1#o(0,1)="}]} :deltas [{:stereo-atom {:modify [0 "Th#o(0,1)/"]}}]}"##)]
    #[case::stereo_atom_ligand_symmetry_removal(r##"{:lhs {:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1#p(0,1)"}]} :deltas [{:stereo-atom {:modify [0 "Th#p(0,1)*"]}}]}"##)]
    fn test_reaction_dsl_from_edn_to_edn_roundtrip(#[case] input: &str) {
        let dsl = ReactionDsl::from_edn(&read_string(input).unwrap()).unwrap();
        let reparsed = ReactionDsl::from_edn(&dsl.to_edn()).unwrap();
        assert_eq!(reparsed, dsl);
    }

    /// Every `fuzz_reaction` seed must parse and satisfy the fuzz invariant (streaming and tree
    /// parsers agree) — guards the seed corpus against rot as the DSL evolves.
    #[rstest]
    fn test_fuzz_reaction_seeds_valid() {
        let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/fuzz/seeds/fuzz_reaction");
        let mut count = 0;
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let data = fs::read_to_string(&path).unwrap();
            let stream = ReactionDsl::from_edn_str(&data).ok();
            let tree = read_string(&data)
                .ok()
                .and_then(|edn| ReactionDsl::from_edn(&edn).ok());
            assert!(stream.is_some(), "seed {path:?} failed to parse");
            assert_eq!(
                stream, tree,
                "seed {path:?}: streaming and tree parsers disagree"
            );
            count += 1;
        }
        assert_eq!(count, 31);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::modify(r##"{:lhs {:atoms ["Br#c0"]} :deltas [{:atom {:modify [0 "#c-1"]}}]}"##)]
    #[case::add_atom_and_bond(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:bond {:add [0 1 "1"]}}]}"##)]
    #[case::molecule_constraint(r##"{:lhs {:atoms ["C" "N"] :bonds [[0 1 "1"]]} :deltas [{:constraint {:add {:connected {}}}}]}"##)]
    #[case::entity_leaf_constraint(r##"{:lhs {:atoms ["C"]} :deltas [{:atom {:add "O"}} {:constraint {:add {:atom [1 {:valence 2}]}}}]}"##)]
    fn test_reaction_ast_to_edn(#[case] input: &str) {
        let ast = ReactionAst::from_edn(&read_string(input).unwrap()).unwrap();
        let reparsed = ReactionAst::from_edn(&ast.to_edn()).unwrap();
        assert_eq!(reparsed, ast);
    }

    /// lhs C(0)–O(1) with bond 0; deltas add N(2) then the bond (1, 2) as delta bond 1.
    #[fixture]
    fn add_bond_reaction() -> ReactionNamespace {
        let input = r##"{:lhs {:atoms ["C" "O"] :bonds [[0 1 "1"]]} :deltas [{:atom {:add [:x "N"]}} {:bond {:add [1 :x "1"]}}]}"##;
        let reaction = ReactionAst::from_edn(&read_string(input).unwrap()).unwrap();
        ReactionNamespace::from_ast(&reaction)
    }

    // `from_ast` reproduces the parse-time namespace across both regions: a structural bond ref to an
    // lhs pair resolves to its lhs id, to a delta-added pair to its delta id, and to a non-pair fails.
    #[rstest]
    #[case::lhs_bond(0, 1, Ok(BondId(0)))]
    #[case::delta_bond(1, 2, Ok(BondId(1)))]
    #[case::self_pair(
        0,
        0,
        Err(ParseError::InvalidRef { kind: "bond", value: "[0 0]".to_string() })
    )]
    fn test_reaction_namespace_from_ast(
        #[from(add_bond_reaction)] ns: ReactionNamespace,
        #[case] a: usize,
        #[case] b: usize,
        #[case] expected: Result<BondId, ParseError>,
    ) {
        let structural = BondRef::Structural([AtomRef::Index(a), AtomRef::Index(b)]);
        assert_eq!(structural.resolve(&ns), expected);
    }
}
