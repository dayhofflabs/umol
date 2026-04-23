//! Molecule DSL.
//!
//! `MoleculeDsl` wraps a `MoleculeAst` together with the `Metadata` that records
//! the surface-form id/alias bindings (atom ids, bond ids, etc.). The EDN
//! form is a map keyed by `:atoms`, `:bonds`, `:dative`, `:aromatic`,
//! `:multicenter`, `:noncovalent`, `:atom-aliases`/`:aliases`, and
//! `:constraints`. Each entity delegates to its own entity DSL. Constraints
//! parse directly into the typed `Constraint` tree.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use bimap::BiMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, FromEdn, ToEdn};

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::BondDsl;
use super::constraint::{
    AromaticSystemRef, AtomRef, BondRef, ConstraintsDsl, DativeBondRef, MulticenterBondRef,
    NoncovalentBondRef, ResolveContext,
};
use super::dative::DativeBondDsl;
use super::error::ParseError;
use super::multicenter::MulticenterBondDsl;
use super::noncovalent::NoncovalentBondDsl;
use crate::ast::atom::AtomAst;
use crate::ast::constraint::{
    AromaticSystemConstraint, AtomConstraint, BondConstraint, DativeBondConstraint,
    MulticenterBondConstraint, NoncovalentBondConstraint,
};
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::molecule::MoleculeAst;
use crate::ast::spin::SpinStateAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;
use crate::dsl::config::MoleculeDefaults;

/// Surface-form metadata paired with a `MoleculeAst`. Records atom ids,
/// per-entity ids, and the atom-alias table. Never drifts: rewrapped
/// atomically through `MoleculeDsl::from_parts`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Metadata {
    pub atom_ids: IndexMap<AtomIdx, String>,
    pub atom_aliases: BiMap<String, AtomDsl>,
    pub bond_ids: IndexMap<BondIdx, String>,
    pub dative_bond_ids: IndexMap<DativeBondIdx, String>,
    pub aromatic_system_ids: IndexMap<AromaticSystemIdx, String>,
    pub multicenter_bond_ids: IndexMap<MulticenterBondIdx, String>,
    pub noncovalent_bond_ids: IndexMap<NoncovalentBondIdx, String>,
}

/// Surface DSL for a whole molecule. Pairs `MoleculeAst` with `Metadata`;
/// fields are private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default)]
pub struct MoleculeDsl {
    ast: MoleculeAst,
    metadata: Metadata,
}

impl MoleculeDsl {
    pub fn from_parts(ast: MoleculeAst, metadata: Metadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.ast
    }

    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (MoleculeAst, Metadata) {
        (self.ast, self.metadata)
    }
}

impl PartialEq for MoleculeDsl {
    fn eq(&self, other: &Self) -> bool {
        self.ast == other.ast && self.metadata == other.metadata
    }
}

impl Eq for MoleculeDsl {}

impl FromStr for MoleculeDsl {
    type Err = ParseError;

    fn from_str(_s: &str) -> Result<Self, Self::Err> {
        todo!("MoleculeDsl::from_str (Phase 5)")
    }
}

impl Display for MoleculeDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn(_edn: &Edn<'de>) -> Result<Self, DeError> {
        todo!("MoleculeDsl::from_edn (Phase 3)")
    }

    fn from_edn_str(_input: &'de str) -> Result<Self, EdnError> {
        todo!("MoleculeDsl::from_edn_str (Phase 4)")
    }
}

impl ToEdn for MoleculeDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(&self.ast, &self.metadata)
    }
}

// -- Render --------------------

fn render_molecule_edn(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(8);
    map.insert(Edn::keyword("atoms"), render_atoms(ast, meta));
    map.insert(Edn::keyword("bonds"), render_bonds(ast, meta));
    if ast.dative_bond_count() > 0 {
        map.insert(Edn::keyword("dative"), render_dative(ast, meta));
    }
    if ast.aromatic_system_count() > 0 {
        map.insert(Edn::keyword("aromatic"), render_aromatic(ast, meta));
    }
    if ast.multicenter_bond_count() > 0 {
        map.insert(Edn::keyword("multicenter"), render_multicenter(ast, meta));
    }
    if ast.noncovalent_bond_count() > 0 {
        map.insert(Edn::keyword("noncovalent"), render_noncovalent(ast, meta));
    }
    if !meta.atom_aliases.is_empty() {
        map.insert(Edn::keyword("atom-aliases"), render_atom_aliases(meta));
    }
    if !ast.constraints().is_empty() {
        let ctx = ResolveContext::for_rendering(meta);
        let dsl = ConstraintsDsl::from_ast(ast.constraints(), &ctx)
            .expect("ConstraintsDsl::from_ast is infallible for a well-formed AST");
        map.insert(Edn::keyword("constraints"), dsl.to_edn());
    }
    Edn::Map(map)
}

fn render_atoms(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .atoms()
        .iter()
        .map(|view| render_atom_entry(view.idx, view.data, meta))
        .collect();
    Edn::Vector(entries.into())
}

fn render_atom_entry(idx: AtomIdx, atom: &AtomAst, meta: &Metadata) -> Edn<'static> {
    let dsl = AtomDsl(atom.clone());
    let spec = if let Some(alias) = meta.atom_aliases.get_by_right(&dsl) {
        Edn::Keyword(EdnKeyword::owned(alias.clone()))
    } else {
        dsl.to_edn()
    };
    match meta.atom_ids.get(&idx) {
        Some(id) => Edn::Vector(
            vec![Edn::Keyword(EdnKeyword::owned(id.clone())), spec].into(),
        ),
        None => spec,
    }
}

fn render_atom_ref(idx: AtomIdx, meta: &Metadata) -> Edn<'static> {
    match meta.atom_ids.get(&idx) {
        Some(id) => Edn::Keyword(EdnKeyword::owned(id.clone())),
        None => Edn::Int(idx.index() as i64),
    }
}

fn render_bonds(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .map(|view| {
            let bond_edn = BondDsl(view.data.clone()).to_edn();
            let a = render_atom_ref(view.src, meta);
            let b = render_atom_ref(view.tgt, meta);
            match meta.bond_ids.get(&view.idx) {
                Some(id) => {
                    let mut m = EdnMap::with_capacity(4);
                    m.insert(
                        Edn::keyword("id"),
                        Edn::Keyword(EdnKeyword::owned(id.clone())),
                    );
                    m.insert(Edn::keyword("a"), a);
                    m.insert(Edn::keyword("b"), b);
                    m.insert(Edn::keyword("type"), bond_edn);
                    Edn::Map(m)
                }
                None => Edn::Vector(vec![a, b, bond_edn].into()),
            }
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_dative(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .dative_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.dative_bond_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            m.insert(Edn::keyword("donor"), render_atom_ref(view.donor, meta));
            m.insert(
                Edn::keyword("acceptor"),
                render_atom_ref(view.acceptor, meta),
            );
            let type_str = DativeBondDsl(view.data.clone()).to_string();
            if !type_str.is_empty() {
                m.insert(Edn::keyword("type"), Edn::Str(Cow::Owned(type_str)));
            }
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_aromatic(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .aromatic_systems()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(3);
            if let Some(id) = meta.aromatic_system_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            let atoms: Vec<Edn<'static>> =
                view.atoms().map(|a| render_atom_ref(a, meta)).collect();
            m.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));
            let type_str = AromaticSystemDsl(view.data.clone()).to_string();
            if !type_str.is_empty() {
                m.insert(Edn::keyword("type"), Edn::Str(Cow::Owned(type_str)));
            }
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_multicenter(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .multicenter_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(3);
            if let Some(id) = meta.multicenter_bond_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            let atoms: Vec<Edn<'static>> =
                view.atoms().map(|a| render_atom_ref(a, meta)).collect();
            m.insert(Edn::keyword("atoms"), Edn::Vector(atoms.into()));
            let type_str = MulticenterBondDsl(view.data.clone()).to_string();
            if !type_str.is_empty() {
                m.insert(Edn::keyword("type"), Edn::Str(Cow::Owned(type_str)));
            }
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_noncovalent(ast: &MoleculeAst, meta: &Metadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .noncovalent_bonds()
        .iter()
        .map(|view| {
            let mut m = EdnMap::with_capacity(4);
            if let Some(id) = meta.noncovalent_bond_ids.get(&view.idx) {
                m.insert(
                    Edn::keyword("id"),
                    Edn::Keyword(EdnKeyword::owned(id.clone())),
                );
            }
            m.insert(Edn::keyword("a"), render_atom_ref(view.atoms[0], meta));
            m.insert(Edn::keyword("b"), render_atom_ref(view.atoms[1], meta));
            m.insert(
                Edn::keyword("type"),
                NoncovalentBondDsl(view.data.clone()).to_edn(),
            );
            Edn::Map(m)
        })
        .collect();
    Edn::Vector(entries.into())
}

fn render_atom_aliases(meta: &Metadata) -> Edn<'static> {
    let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(meta.atom_aliases.len() * 2);
    for (name, dsl) in meta.atom_aliases.iter() {
        pairs.push(Edn::Keyword(EdnKeyword::owned(name.clone())));
        pairs.push(dsl.to_edn());
    }
    Edn::Vector(pairs.into())
}

impl FromAst<MoleculeAst> for MoleculeDsl {
    type Ctx<'a> = MoleculeDefaults;
    type Error = ParseError;

    fn from_ast<'a>(
        _ast: &MoleculeAst,
        _cfg: &Self::Ctx<'a>,
    ) -> Result<Self, ParseError> {
        todo!("MoleculeDsl::from_ast (Phase 5)")
    }
}

impl IntoAst<MoleculeAst> for MoleculeDsl {
    type Ctx<'a> = MoleculeDefaults;
    type Error = ParseError;

    fn into_ast<'a>(self, _cfg: &Self::Ctx<'a>) -> Result<MoleculeAst, ParseError> {
        todo!("MoleculeDsl::to_ast (Phase 5)")
    }
}

// -- Private parse intermediate ---------------------------------------------
//
// Unresolved, owned-by-value tree that mirrors the EDN shape. Atom entries and
// per-bond endpoints carry `AtomRef` (index or id); constraint leaves carry
// typed per-entity `Constraint*` variants already parsed from their EDN form.
// Lowered destructively via `into_ast(self, cfg)` so that allocations move
// into the final `MoleculeAst`.

#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct MoleculeInput {
    pub(crate) atoms: Vec<AtomEntryInput>,
    pub(crate) bonds: Vec<BondEntryInput>,
    pub(crate) dative_bonds: Vec<DativeBondEntryInput>,
    pub(crate) aromatic_systems: Vec<AromaticSystemEntryInput>,
    pub(crate) multicenter_bonds: Vec<MulticenterBondEntryInput>,
    pub(crate) noncovalent_bonds: Vec<NoncovalentBondEntryInput>,
    pub(crate) atom_aliases: Vec<(String, AtomDsl)>,
    pub(crate) constraints: Vec<ConstraintInput>,
}

/// Atom entry in a parsed molecule map. Mirrors the DSL spec §4 grammar
/// `atom-entry ::= atom-spec | [ keyword atom-spec ]`.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AtomEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) spec: AtomSpecInput,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum AtomSpecInput {
    Bare(Box<AtomDsl>),
    Alias(String),
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct BondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct DativeBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) donor: AtomRef,
    pub(crate) acceptor: AtomRef,
    pub(crate) bond: DativeBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct AromaticSystemEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) system: AromaticSystemDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct MulticenterBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) bond: MulticenterBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct NoncovalentBondEntryInput {
    pub(crate) id: Option<String>,
    pub(crate) a: AtomRef,
    pub(crate) b: AtomRef,
    pub(crate) bond: NoncovalentBondDsl,
}

/// Pre-resolution mirror of `Constraint`. Entity leaves carry the unresolved
/// reference; combinators recurse; molecule-scope leaves are parsed in place.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum ConstraintInput {
    Atom(AtomRef, AtomConstraint),
    Bond(BondRef, BondConstraint),
    DativeBond(DativeBondRef, DativeBondConstraint),
    AromaticSystem(AromaticSystemRef, AromaticSystemConstraint),
    MulticenterBond(MulticenterBondRef, MulticenterBondConstraint),
    NoncovalentBond(NoncovalentBondRef, NoncovalentBondConstraint),
    Molecule(MoleculeConstraintInput),
    And(Vec<ConstraintInput>),
    Or(Vec<ConstraintInput>),
    Not(Box<ConstraintInput>),
}

/// Pre-resolution mirror of `MoleculeConstraint`. `SubPattern.pattern` holds a
/// nested `MoleculeInput` so the sub-molecule is resolved against its own
/// local id scope before the outer molecule completes lowering.
#[derive(Clone, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum MoleculeConstraintInput {
    ChargeSum {
        atoms: Vec<AtomRef>,
        sum: ValueAst,
    },
    SpinSum {
        atoms: Vec<AtomRef>,
        spin: SpinStateAst,
    },
    BondOrderSum {
        bonds: Vec<BondRef>,
        sum: ValueAst,
    },
    Connected(Vec<AtomRef>),
    SubPattern {
        anchor: SubPatternAnchorInput,
        pattern: Box<MoleculeInput>,
    },
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) struct SubPatternAnchorInput {
    pub(crate) atoms: Vec<(AtomRef, AtomRef)>,
    pub(crate) bonds: Vec<(BondRef, BondRef)>,
    pub(crate) dative_bonds: Vec<(DativeBondRef, DativeBondRef)>,
    pub(crate) aromatic_systems: Vec<(AromaticSystemRef, AromaticSystemRef)>,
    pub(crate) multicenter_bonds: Vec<(MulticenterBondRef, MulticenterBondRef)>,
    pub(crate) noncovalent_bonds: Vec<(NoncovalentBondRef, NoncovalentBondRef)>,
}

impl MoleculeInput {
    /// Destructive lowering: consumes the input, resolves refs against the
    /// built id scopes, lifts bare entity-leaf constraints onto their entity
    /// AST's inline store, and produces the final `MoleculeAst` with its
    /// `Metadata`. Called from `FromEdn::from_edn` and the streaming path.
    #[allow(dead_code)]
    pub(crate) fn into_ast(
        self,
        _cfg: &MoleculeDefaults,
    ) -> Result<(MoleculeAst, Metadata), ParseError> {
        todo!("MoleculeInput::into_ast (Phase 3)")
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::atom::{AtomAst, ElementAst};
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::Constraints;

    fn c_atom() -> AtomAst {
        AtomAst::new(ElementAst::Lit(Element::C))
    }

    fn single_bond() -> BondAst {
        BondAst::new(ValueAst::Lit(1))
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_empty() {
        let ast = MoleculeAst::default();
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        assert_eq!(edn, read_string("{:atoms [] :bonds []}").unwrap());
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_two_atoms_one_bond() {
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_atom_with_id() {
        let mut atom_ids = IndexMap::new();
        atom_ids.insert(AtomIdx(0), "c1".to_string());
        let meta = Metadata {
            atom_ids,
            ..Metadata::default()
        };
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, meta);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms [[:c1 "C"] "C"] :bonds []}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_bond_with_id_uses_map_form() {
        let mut bond_ids = IndexMap::new();
        bond_ids.insert(BondIdx(0), "b1".to_string());
        let meta = Metadata {
            bond_ids,
            ..Metadata::default()
        };
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, meta);
        let edn = dsl.to_edn();
        assert_eq!(
            edn,
            read_string(r##"{:atoms ["C" "C"] :bonds [{:id :b1 :a 0 :b 1 :type "1"}]}"##)
                .unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_atom_alias_substituted() {
        let mut atom_aliases: BiMap<String, AtomDsl> = BiMap::new();
        atom_aliases.insert("x".into(), AtomDsl(c_atom()));
        let meta = Metadata {
            atom_aliases,
            ..Metadata::default()
        };
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, meta);
        let edn = dsl.to_edn();
        // Both atoms match the alias — rendered as :x keyword references; the
        // alias table emits the :atom-aliases key.
        assert_eq!(
            edn,
            read_string(r##"{:atoms [:x :x] :bonds [] :atom-aliases [:x "C"]}"##).unwrap()
        );
    }

    #[rstest]
    fn test_molecule_dsl_display_matches_edn() {
        let ast = MoleculeAst::new(
            vec![c_atom(), c_atom()],
            vec![(AtomIdx(0), AtomIdx(1), single_bond())],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        assert_eq!(dsl.to_string(), dsl.to_edn().to_string());
    }

    #[rstest]
    fn test_molecule_dsl_to_edn_omits_empty_optional_sections() {
        let ast = MoleculeAst::new(
            vec![c_atom()],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::new(),
        );
        let dsl = MoleculeDsl::from_parts(ast, Metadata::default());
        let edn = dsl.to_edn();
        let Edn::Map(m) = &edn else {
            panic!("expected map");
        };
        assert!(m.get_keyword("dative").is_none());
        assert!(m.get_keyword("aromatic").is_none());
        assert!(m.get_keyword("multicenter").is_none());
        assert!(m.get_keyword("noncovalent").is_none());
        assert!(m.get_keyword("atom-aliases").is_none());
        assert!(m.get_keyword("constraints").is_none());
    }
}
