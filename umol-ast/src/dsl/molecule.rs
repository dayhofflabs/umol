//! Molecule DSL.
//!
//! `MoleculeDsl` wraps a `MoleculeAst` together with the `MoleculeMetadata` that records
//! the surface-form entity keywords and atom aliases. The EDN
//! form is a map keyed by `:atoms`, `:bonds`, `:dative-bonds`, `:aromatic-systems`,
//! `:multicenter-bonds`, `:noncovalent-bonds`, `:atom-aliases`/`:aliases`, and
//! `:constraints`. Each entity delegates to its own entity DSL. Constraints
//! parse directly into the typed `Constraint` tree.

// Closures like `|e| T::from_edn(e)` passed to `parse_vec` can't be replaced
// by bare `T::from_edn` — type-erasing the fn item loses the `for<'a>` HRTB
// on the `FromEdn<'a>` impl.
#![allow(clippy::redundant_closure)]

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};

use super::aromatic::AromaticSystemDsl;
use super::atom::AtomDsl;
use super::bond::{expand_bond_keyword, BondDsl};
use super::config::MoleculeDefaults;
use super::constraint::{read_constraints_dsl, ConstraintDsl, ConstraintsDsl};
use super::dative::DativeBondDsl;
use super::edn_utils::{
    atoms_pair, atoms_vec, eof_err, missing, optional_id_keyword, parse_vec, read_map, read_vec,
    required_key, two_atom_refs, unexpected_byte_kind,
};
use super::error::ParseError;
use super::metadata::{Metadata, MetadataError, MoleculeMetadata};
use super::multicenter::MulticenterBondDsl;
use super::namespace::{MoleculeContext, Namespace};
use super::noncovalent::NoncovalentBondDsl;
use super::refs::{
    parse_stereo_ligand, read_atom_ref, read_bond_ref, read_stereo_ligand, AtomRef, BondRef,
    StereoLigandRef,
};
use super::stereo::{
    expand_stereo_atom_keyword, expand_stereo_bond_keyword, StereoAtomDsl, StereoBondDsl,
};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::dative::DativeBondAst;
use crate::ast::entity::Entity;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ast::ligand::{StereoLigand, StereoLigandKind};
use crate::ast::molecule::{MoleculeAst, MoleculeParts};
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::noncovalent::NoncovalentBondAst;
use crate::ast::stereo::{StereoAtomAst, StereoBondAst};
use crate::ast::traits::{FromAst, IntoAst};

/// Surface DSL for a whole molecule. Pairs `MoleculeAst` with `MoleculeMetadata`;
/// fields are private so metadata cannot drift onto a different AST.
#[derive(Clone, Debug, Default)]
pub struct MoleculeDsl {
    ast: MoleculeAst,
    metadata: MoleculeMetadata,
}

impl MoleculeDsl {
    /// Pair a molecule AST with coherent surface metadata.
    pub fn new(ast: MoleculeAst, metadata: MoleculeMetadata) -> Result<Self, MetadataError> {
        for (entity, _) in metadata.iter_keywords() {
            let contains = match entity {
                Entity::Atom(id) => ast.atoms().contains(id),
                Entity::Bond(id) => ast.bonds().contains(id),
                Entity::DativeBond(id) => ast.dative_bonds().contains(id),
                Entity::AromaticSystem(id) => ast.aromatic_systems().contains(id),
                Entity::MulticenterBond(id) => ast.multicenter_bonds().contains(id),
                Entity::NoncovalentBond(id) => ast.noncovalent_bonds().contains(id),
                Entity::StereoAtom(id) => ast.stereo_atoms().contains(id),
                Entity::StereoBond(id) => ast.stereo_bonds().contains(id),
            };
            if !contains {
                return Err(MetadataError::EntityOutOfRange(entity));
            }
        }
        Ok(Self::from_parts(ast, metadata))
    }

    fn from_parts(ast: MoleculeAst, metadata: MoleculeMetadata) -> Self {
        Self { ast, metadata }
    }

    pub fn ast(&self) -> &MoleculeAst {
        &self.ast
    }

    pub fn metadata(&self) -> &MoleculeMetadata {
        &self.metadata
    }

    pub fn into_parts(self) -> (MoleculeAst, MoleculeMetadata) {
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

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        MoleculeDsl::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for MoleculeDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let input = parse_molecule_input(edn)?;
        let (ast, context) = input
            .into_ast()
            .map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(ast, context.into_metadata()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let mut de = EdnStreamDeserializer::new(input);
        let mi = read_molecule_input(&mut de)?;
        de.expect_eof()?;
        let (ast, context) = mi.into_ast().map_err(|e| DeError::Custom(e.to_string()))?;
        Ok(MoleculeDsl::from_parts(ast, context.into_metadata()))
    }
}

/// Direct EDN parsing for `MoleculeAst`. Accepts the same molecule-map
/// surface as [`MoleculeDsl::from_edn`]; any entity keywords or aliases in the
/// input resolve to positional indices, then the metadata is discarded —
/// the result is metadata-free.
impl<'de> FromEdn<'de> for MoleculeAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        MoleculeDsl::from_edn(edn).map(|dsl| dsl.into_parts().0)
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        MoleculeDsl::from_edn_str(input).map(|dsl| dsl.into_parts().0)
    }
}

impl FromStr for MoleculeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_edn_str(s).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

impl Display for MoleculeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_edn())
    }
}

/// Direct EDN rendering for `MoleculeAst`. Always emits canonical positional
/// refs (no entity keywords, no aliases) since the AST carries no metadata.
/// For keyword-bearing surface output, wrap in [`MoleculeDsl`] with appropriate
/// [`MoleculeMetadata`] and call [`MoleculeDsl::to_edn`].
impl ToEdn for MoleculeAst {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(self, &MoleculeMetadata::default())
    }
}

// Streaming parse of the molecule map.
pub(super) fn read_molecule_input(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MoleculeInput, EdnError> {
    de.consume_byte(b'{')?;
    let mut mi = MoleculeInput::default();
    loop {
        if de.try_consume_byte(b'}')? {
            break;
        }
        let key = de.read_keyword_name()?;
        match key.as_ref() {
            "atoms" => mi.atoms = read_vec(de, read_atom_entry)?,
            "bonds" => mi.bonds = read_vec(de, read_bond_entry)?,
            "dative-bonds" => mi.dative_bonds = read_vec(de, read_dative_bond_entry)?,
            "aromatic-systems" => mi.aromatic_systems = read_vec(de, read_aromatic_system_entry)?,
            "multicenter-bonds" => {
                mi.multicenter_bonds = read_vec(de, read_multicenter_bond_entry)?
            }
            "noncovalent-bonds" => {
                mi.noncovalent_bonds = read_vec(de, read_noncovalent_bond_entry)?
            }
            "stereo-atoms" => mi.stereo_atoms = read_vec(de, read_stereo_atom_entry)?,
            "stereo-bonds" => mi.stereo_bonds = read_vec(de, read_stereo_bond_entry)?,
            "atom-aliases" => mi.atom_aliases = read_atom_aliases(de)?,
            "constraints" => mi.constraints = read_constraints_dsl(de)?,
            "guards" => {
                de.read_skip_value()?;
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["molecule".into()],
                }
                .into());
            }
        }
    }
    Ok(mi)
}

pub(super) fn read_atom_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AtomEntryInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'[' => {
            de.consume_byte(b'[')?;
            let keyword = de.read_keyword_name()?.into_owned();
            let spec = read_atom_spec(de)?;
            de.consume_byte(b']')?;
            Ok(AtomEntryInput {
                keyword: Some(keyword),
                spec,
            })
        }
        _ => Ok(AtomEntryInput {
            keyword: None,
            spec: read_atom_spec(de)?,
        }),
    }
}

fn read_atom_spec(de: &mut EdnStreamDeserializer<'_>) -> Result<AtomSpecInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'"' => {
            let s = de.read_string()?;
            let dsl: AtomDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("atom", e))?;
            Ok(AtomSpecInput::Bare(Box::new(dsl)))
        }
        b':' => {
            let name = de.read_keyword_name()?;
            Ok(AtomSpecInput::Alias(name.into_owned()))
        }
        b => Err(DeError::TypeMismatch {
            expected: "atom-string or :alias",
            got: unexpected_byte_kind(b),
            path: vec!["atom-spec".into()],
        }
        .into()),
    }
}

fn read_bond_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<BondDsl, EdnError> {
    let byte = de.peek_byte()?.ok_or_else(eof_err)?;
    let text: Cow<'_, str> = match byte {
        b':' => {
            let name = de.read_keyword_name()?;
            let expanded = expand_bond_keyword(name.as_ref())
                .ok_or_else(|| DeError::Custom(format!("unknown bond keyword :{}", name)))?;
            Cow::Borrowed(expanded)
        }
        _ => de.read_string()?,
    };
    text.as_ref()
        .parse()
        .map_err(|e| DeError::subgrammar("bond", e).into())
}

fn read_dative_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<DativeBondDsl, EdnError> {
    let byte = de.peek_byte()?.ok_or_else(eof_err)?;
    let text: Cow<'_, str> = match byte {
        b':' => {
            let name = de.read_keyword_name()?;
            let expanded = super::dative::expand_dative_keyword(name.as_ref())
                .ok_or_else(|| DeError::Custom(format!("unknown dative keyword :{}", name)))?;
            Cow::Borrowed(expanded)
        }
        _ => de.read_string()?,
    };
    text.as_ref()
        .parse()
        .map_err(|e| DeError::subgrammar("dative", e).into())
}

/// A two-endpoint `:atoms` vector for a binary relation: exactly two refs.
pub(super) fn read_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<BondEntryInput, EdnError> {
    match de.peek_byte()?.ok_or_else(eof_err)? {
        b'[' => {
            de.consume_byte(b'[')?;
            let a = read_atom_ref(de)?;
            let b = read_atom_ref(de)?;
            let bond = read_bond_dsl(de)?;
            de.consume_byte(b']')?;
            Ok(BondEntryInput {
                keyword: None,
                first: a,
                second: b,
                bond,
            })
        }
        b'{' => {
            let mut keyword = None;
            let mut atoms = None;
            let mut bond = None;
            read_map(de, |de, key| {
                match key {
                    "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
                    "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
                    "type" => bond = Some(read_bond_dsl(de)?),
                    _ => de.read_skip_value()?,
                }
                Ok(())
            })?;
            let atoms = atoms.ok_or_else(|| missing("atoms", "bond-entry"))?;
            let [a, b] = two_atom_refs(atoms, "bond-entry")?;
            Ok(BondEntryInput {
                keyword,
                first: a,
                second: b,
                bond: bond.ok_or_else(|| missing("type", "bond-entry"))?,
            })
        }
        bb => Err(DeError::TypeMismatch {
            expected: "bond-entry map or 3-vec",
            got: unexpected_byte_kind(bb),
            path: vec!["bond-entry".into()],
        }
        .into()),
    }
}

pub(super) fn read_dative_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<DativeBondEntryInput, EdnError> {
    let mut keyword = None;
    let mut donors = None;
    let mut acceptor = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
            "donors" => donors = Some(read_vec(de, read_atom_ref)?),
            "acceptor" => acceptor = Some(read_atom_ref(de)?),
            "type" => bond = Some(read_dative_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(DativeBondEntryInput {
        keyword,
        donors: donors.ok_or_else(|| missing("donors", "dative-bond-entry"))?,
        acceptor: acceptor.ok_or_else(|| missing("acceptor", "dative-bond-entry"))?,
        bond: bond.ok_or_else(|| missing("type", "dative-bond-entry"))?,
    })
}

pub(super) fn read_aromatic_system_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<AromaticSystemEntryInput, EdnError> {
    let mut keyword = None;
    let mut atoms = None;
    let mut system = None;
    read_map(de, |de, key| {
        match key {
            "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "type" => {
                let s = de.read_string()?;
                system = Some(
                    s.as_ref()
                        .parse::<AromaticSystemDsl>()
                        .map_err(|e| DeError::subgrammar("aromatic", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    let system = system.ok_or_else(|| missing("type", "aromatic-system-entry"))?;
    Ok(AromaticSystemEntryInput {
        keyword,
        atoms: atoms.ok_or_else(|| missing("atoms", "aromatic-system-entry"))?,
        system,
    })
}

pub(super) fn read_multicenter_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<MulticenterBondEntryInput, EdnError> {
    let mut keyword = None;
    let mut atoms = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "type" => {
                let s = de.read_string()?;
                bond = Some(
                    s.as_ref()
                        .parse::<MulticenterBondDsl>()
                        .map_err(|e| DeError::subgrammar("multicenter", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    let bond = bond.ok_or_else(|| missing("type", "multicenter-bond-entry"))?;
    Ok(MulticenterBondEntryInput {
        keyword,
        atoms: atoms.ok_or_else(|| missing("atoms", "multicenter-bond-entry"))?,
        bond,
    })
}

pub(super) fn read_noncovalent_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<NoncovalentBondEntryInput, EdnError> {
    let mut keyword = None;
    let mut atoms = None;
    let mut bond = None;
    read_map(de, |de, key| {
        match key {
            "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
            "atoms" => atoms = Some(read_vec(de, read_atom_ref)?),
            "type" => {
                let text = de.read_string_or_keyword()?;
                bond = Some(
                    text.as_ref()
                        .parse::<NoncovalentBondDsl>()
                        .map_err(|e| DeError::subgrammar("noncovalent", e))?,
                );
            }
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    let atoms = atoms.ok_or_else(|| missing("atoms", "noncovalent-bond-entry"))?;
    let [a, b] = two_atom_refs(atoms, "noncovalent-bond-entry")?;
    Ok(NoncovalentBondEntryInput {
        keyword,
        first: a,
        second: b,
        bond: bond.ok_or_else(|| missing("type", "noncovalent-bond-entry"))?,
    })
}

fn read_stereo_atom_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<StereoAtomDsl, EdnError> {
    if de.peek_byte()?.ok_or_else(eof_err)? == b':' {
        let kw = de.read_keyword_name()?;
        let expanded = expand_stereo_atom_keyword(kw.as_ref())
            .ok_or_else(|| DeError::Custom(format!("unknown stereo atom keyword :{kw}")))?;
        expanded
            .parse::<StereoAtomDsl>()
            .map_err(|e| DeError::subgrammar("stereo atom", e).into())
    } else {
        de.read_string()?
            .as_ref()
            .parse::<StereoAtomDsl>()
            .map_err(|e| DeError::subgrammar("stereo atom", e).into())
    }
}

fn read_stereo_bond_dsl(de: &mut EdnStreamDeserializer<'_>) -> Result<StereoBondDsl, EdnError> {
    if de.peek_byte()?.ok_or_else(eof_err)? == b':' {
        let kw = de.read_keyword_name()?;
        let expanded = expand_stereo_bond_keyword(kw.as_ref())
            .ok_or_else(|| DeError::Custom(format!("unknown stereo bond keyword :{kw}")))?;
        expanded
            .parse::<StereoBondDsl>()
            .map_err(|e| DeError::subgrammar("stereo bond", e).into())
    } else {
        de.read_string()?
            .as_ref()
            .parse::<StereoBondDsl>()
            .map_err(|e| DeError::subgrammar("stereo bond", e).into())
    }
}

pub(super) fn read_stereo_atom_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoAtomEntryInput, EdnError> {
    let mut keyword = None;
    let mut site = None;
    let mut ligands = None;
    let mut stereo = None;
    read_map(de, |de, key| {
        match key {
            "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
            "site" => site = Some(read_atom_ref(de)?),
            "ligands" => ligands = Some(read_vec(de, read_stereo_ligand)?),
            "type" => stereo = Some(read_stereo_atom_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(StereoAtomEntryInput {
        keyword,
        site: site.ok_or_else(|| missing("site", "stereo-atom-entry"))?,
        ligands: ligands.ok_or_else(|| missing("ligands", "stereo-atom-entry"))?,
        stereo: stereo.ok_or_else(|| missing("type", "stereo-atom-entry"))?,
    })
}

pub(super) fn read_stereo_bond_entry(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<StereoBondEntryInput, EdnError> {
    let mut keyword = None;
    let mut site = None;
    let mut ligands = None;
    let mut stereo = None;
    read_map(de, |de, key| {
        match key {
            "id" => keyword = Some(de.read_keyword_name()?.into_owned()),
            "site" => site = Some(read_bond_ref(de)?),
            "ligands" => ligands = Some(read_vec(de, read_stereo_ligand)?),
            "type" => stereo = Some(read_stereo_bond_dsl(de)?),
            _ => de.read_skip_value()?,
        }
        Ok(())
    })?;
    Ok(StereoBondEntryInput {
        keyword,
        site: site.ok_or_else(|| missing("site", "stereo-bond-entry"))?,
        ligands: ligands.ok_or_else(|| missing("ligands", "stereo-bond-entry"))?,
        stereo: stereo.ok_or_else(|| missing("type", "stereo-bond-entry"))?,
    })
}

pub(super) fn read_atom_aliases(
    de: &mut EdnStreamDeserializer<'_>,
) -> Result<Vec<(String, Box<AtomDsl>)>, EdnError> {
    de.consume_byte(b'[')?;
    let mut out = Vec::new();
    loop {
        if de.try_consume_byte(b']')? {
            break;
        }
        let name = de.read_keyword_name()?.into_owned();
        if de.try_consume_byte(b']')? {
            return Err(DeError::Custom(
                ":atom-aliases must have even length (keyword/atom-string pairs)".into(),
            )
            .into());
        }
        let s = de.read_string()?;
        let dsl: AtomDsl = s
            .as_ref()
            .parse()
            .map_err(|e| DeError::subgrammar("atom", e))?;
        out.push((name, Box::new(dsl)));
    }
    Ok(out)
}

impl ToEdn for MoleculeDsl {
    fn to_edn(&self) -> Edn<'static> {
        render_molecule_edn(&self.ast, &self.metadata)
    }
}

impl FromAst<MoleculeAst> for MoleculeDsl {
    type Ctx = MoleculeDefaults;

    fn from_ast(ast: &MoleculeAst, cfg: &Self::Ctx) -> Self {
        let mut ast_out = ast.clone();
        ast_out.modify_atoms(|atom| AtomDsl::from_ast(&atom, &cfg.atom).0);
        ast_out.modify_bonds(|bond| BondDsl::from_ast(&bond, &cfg.bond).0);
        ast_out.modify_aromatic_systems(|system| {
            AromaticSystemDsl::from_ast(&system, &cfg.aromatic_system).0
        });
        ast_out.modify_multicenter_bonds(|bond| {
            MulticenterBondDsl::from_ast(&bond, &cfg.multicenter_bond).0
        });
        ast_out.modify_dative_bonds(|bond| DativeBondDsl::from_ast(&bond, &cfg.dative_bond).0);
        ast_out.modify_noncovalent_bonds(|bond| {
            NoncovalentBondDsl::from_ast(&bond, &cfg.noncovalent_bond).0
        });
        ast_out.modify_stereo_atoms(|stereo_atom| {
            StereoAtomDsl::from_ast(&stereo_atom, &cfg.stereo_atom).0
        });
        ast_out.modify_stereo_bonds(|stereo_bond| {
            StereoBondDsl::from_ast(&stereo_bond, &cfg.stereo_bond).0
        });
        MoleculeDsl {
            ast: ast_out,
            metadata: MoleculeMetadata::default(),
        }
    }
}

impl IntoAst<MoleculeAst> for MoleculeDsl {
    type Ctx = MoleculeDefaults;

    fn into_ast(self, cfg: &Self::Ctx) -> MoleculeAst {
        let mut ast = self.ast;
        ast.modify_atoms(|atom| AtomDsl(atom).into_ast(&cfg.atom));
        ast.modify_bonds(|bond| BondDsl(bond).into_ast(&cfg.bond));
        ast.modify_dative_bonds(|bond| DativeBondDsl(bond).into_ast(&cfg.dative_bond));
        ast.modify_aromatic_systems(|system| {
            AromaticSystemDsl(system).into_ast(&cfg.aromatic_system)
        });
        ast.modify_multicenter_bonds(|bond| {
            MulticenterBondDsl(bond).into_ast(&cfg.multicenter_bond)
        });
        ast.modify_noncovalent_bonds(|bond| {
            NoncovalentBondDsl(bond).into_ast(&cfg.noncovalent_bond)
        });
        ast.modify_stereo_atoms(|stereo_atom| {
            StereoAtomDsl(stereo_atom).into_ast(&cfg.stereo_atom)
        });
        ast.modify_stereo_bonds(|stereo_bond| {
            StereoBondDsl(stereo_bond).into_ast(&cfg.stereo_bond)
        });
        ast
    }
}

pub(super) fn render_molecule_edn(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let mut map = EdnMap::with_capacity(8);
    map.insert(Edn::keyword("atoms"), render_atoms(ast, meta));
    map.insert(Edn::keyword("bonds"), render_bonds(ast, meta));
    if ast.dative_bonds().count() > 0 {
        map.insert(Edn::keyword("dative-bonds"), render_dative(ast, meta));
    }
    if ast.aromatic_systems().count() > 0 {
        map.insert(Edn::keyword("aromatic-systems"), render_aromatic(ast, meta));
    }
    if ast.multicenter_bonds().count() > 0 {
        map.insert(
            Edn::keyword("multicenter-bonds"),
            render_multicenter(ast, meta),
        );
    }
    if ast.noncovalent_bonds().count() > 0 {
        map.insert(
            Edn::keyword("noncovalent-bonds"),
            render_noncovalent(ast, meta),
        );
    }
    if ast.stereo_atoms().count() > 0 {
        map.insert(Edn::keyword("stereo-atoms"), render_stereo_atoms(ast, meta));
    }
    if ast.stereo_bonds().count() > 0 {
        map.insert(Edn::keyword("stereo-bonds"), render_stereo_bonds(ast, meta));
    }
    if meta.iter_atom_aliases().len() != 0 {
        map.insert(Edn::keyword("atom-aliases"), render_atom_aliases(meta));
    }
    let constraints_dsl = ConstraintsDsl::from_ast(ast.constraints(), meta)
        .expect("ConstraintsDsl::from_ast is infallible for a well-formed AST");
    if !constraints_dsl.0.is_empty() {
        map.insert(Edn::keyword("constraints"), constraints_dsl.to_edn());
    }
    Edn::Map(map)
}

fn render_atoms(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .atoms()
        .iter()
        .map(|view| render_atom_entry(view.id, view.ast, meta))
        .collect();
    Edn::Vector(entries.into())
}

/// An atom value: its alias keyword if one is bound, else the atom-string.
pub(super) fn render_atom_value(atom: &AtomAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let dsl = AtomDsl::from_ref(atom);
    match meta.atom_alias_name(dsl) {
        Some(alias) => Edn::Keyword(EdnKeyword::owned(alias.to_string())),
        None => dsl.to_edn(),
    }
}

fn render_atom_entry(id: AtomId, atom: &AtomAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let spec = render_atom_value(atom, meta);
    match meta.keyword(Entity::Atom(id)) {
        Some(keyword) => {
            Edn::Vector(vec![Edn::Keyword(EdnKeyword::owned(keyword.to_string())), spec].into())
        }
        None => spec,
    }
}

fn render_atom_ref(id: AtomId, meta: &impl Metadata) -> Edn<'static> {
    match meta.keyword(Entity::Atom(id)) {
        Some(keyword) => Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

/// A bond entry — `[a b type]`, or `{:id … :atoms [a b] :type type}` when the bond has a keyword.
/// `type_edn` is the already-rendered `:type` Edn — one bond-dsl for a molecule, or a `[left right]`
/// vector / op-wrapped map for a span entry; it is not an ast.
pub(super) fn render_bond_entry(
    id: BondId,
    [a, b]: [AtomId; 2],
    type_edn: Edn<'static>,
    meta: &MoleculeMetadata,
) -> Edn<'static> {
    let first = render_atom_ref(a, meta);
    let second = render_atom_ref(b, meta);
    match meta.keyword(Entity::Bond(id)) {
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
            m.insert(Edn::keyword("type"), type_edn);
            Edn::Map(m)
        }
        None => Edn::Vector(vec![first, second, type_edn].into()),
    }
}

fn render_bonds(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .bonds()
        .iter()
        .map(|view| {
            render_bond_entry(
                view.id,
                view.atom_ids(),
                BondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

// Overlay entries: `render_<entity>_entry` builds one entry map (`:id`? + participants + `:type`),
// with `:type` = the caller-supplied `type_edn`. Molecule and reaction-delta rendering pass an
// entity DSL value; reaction-span rendering passes a `{:add|:modify|:remove}`-wrapped value.

pub(super) fn render_dative_entry(
    id: DativeBondId,
    donors: impl Iterator<Item = AtomId>,
    acceptor: AtomId,
    type_edn: Edn<'static>,
    meta: &impl Metadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(keyword) = meta.keyword(Entity::DativeBond(id)) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        );
    }
    m.insert(
        Edn::keyword("donors"),
        Edn::Vector(
            donors
                .map(|a| render_atom_ref(a, meta))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(Edn::keyword("acceptor"), render_atom_ref(acceptor, meta));
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_dative(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .dative_bonds()
        .iter()
        .map(|view| {
            render_dative_entry(
                view.id,
                view.donor_ids(),
                view.acceptor_id(),
                DativeBondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_aromatic_entry(
    id: AromaticSystemId,
    atoms: impl Iterator<Item = AtomId>,
    type_edn: Edn<'static>,
    meta: &impl Metadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    if let Some(keyword) = meta.keyword(Entity::AromaticSystem(id)) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        );
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(
            atoms
                .map(|a| render_atom_ref(a, meta))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_aromatic(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .aromatic_systems()
        .iter()
        .map(|view| {
            render_aromatic_entry(
                view.id,
                view.atom_ids(),
                Edn::Str(Cow::Owned(
                    AromaticSystemDsl::from_ref(view.ast).to_string(),
                )),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_multicenter_entry(
    id: MulticenterBondId,
    atoms: impl Iterator<Item = AtomId>,
    type_edn: Edn<'static>,
    meta: &impl Metadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    if let Some(keyword) = meta.keyword(Entity::MulticenterBond(id)) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        );
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(
            atoms
                .map(|a| render_atom_ref(a, meta))
                .collect::<Vec<_>>()
                .into(),
        ),
    );
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_multicenter(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .multicenter_bonds()
        .iter()
        .map(|view| {
            render_multicenter_entry(
                view.id,
                view.atom_ids(),
                Edn::Str(Cow::Owned(
                    MulticenterBondDsl::from_ref(view.ast).to_string(),
                )),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_noncovalent_entry(
    id: NoncovalentBondId,
    [a, b]: [AtomId; 2],
    type_edn: Edn<'static>,
    meta: &impl Metadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    if let Some(keyword) = meta.keyword(Entity::NoncovalentBond(id)) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        );
    }
    m.insert(
        Edn::keyword("atoms"),
        Edn::Vector(vec![render_atom_ref(a, meta), render_atom_ref(b, meta)].into()),
    );
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_noncovalent(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .noncovalent_bonds()
        .iter()
        .map(|view| {
            render_noncovalent_entry(
                view.id,
                view.atom_ids(),
                NoncovalentBondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_stereo_atom_entry(
    id: StereoAtomId,
    site: AtomId,
    ligands: Vec<Edn<'static>>,
    type_edn: Edn<'static>,
    meta: &impl Metadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(keyword) = meta.keyword(Entity::StereoAtom(id)) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        );
    }
    m.insert(Edn::keyword("site"), render_atom_ref(site, meta));
    m.insert(Edn::keyword("ligands"), Edn::Vector(ligands.into()));
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_stereo_atoms(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .stereo_atoms()
        .iter()
        .map(|view| {
            render_stereo_atom_entry(
                view.id,
                view.site_id(),
                view.ligand_frame()
                    .into_iter()
                    .map(|l| render_stereo_ligand(l, meta))
                    .collect(),
                StereoAtomDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_stereo_bond_entry(
    id: StereoBondId,
    site: BondId,
    ligands: Vec<Edn<'static>>,
    type_edn: Edn<'static>,
    meta: &impl Metadata,
) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(4);
    if let Some(keyword) = meta.keyword(Entity::StereoBond(id)) {
        m.insert(
            Edn::keyword("id"),
            Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        );
    }
    m.insert(Edn::keyword("site"), render_bond_ref(site, meta));
    m.insert(Edn::keyword("ligands"), Edn::Vector(ligands.into()));
    m.insert(Edn::keyword("type"), type_edn);
    Edn::Map(m)
}

fn render_stereo_bonds(ast: &MoleculeAst, meta: &MoleculeMetadata) -> Edn<'static> {
    let entries: Vec<Edn<'static>> = ast
        .stereo_bonds()
        .iter()
        .map(|view| {
            render_stereo_bond_entry(
                view.id,
                view.site_id(),
                view.ligand_frame()
                    .into_iter()
                    .map(|l| render_stereo_ligand(l, meta))
                    .collect(),
                StereoBondDsl::from_ref(view.ast).to_edn(),
                meta,
            )
        })
        .collect();
    Edn::Vector(entries.into())
}

pub(super) fn render_stereo_ligand(ligand: StereoLigand, meta: &impl Metadata) -> Edn<'static> {
    let atom = render_atom_ref(ligand.atom_id, meta);
    match ligand.kind {
        StereoLigandKind::Atom => atom,
        StereoLigandKind::ImplicitHydrogen => Edn::Vector(vec![Edn::keyword("h"), atom].into()),
        StereoLigandKind::LonePair => Edn::Vector(vec![Edn::keyword("lp"), atom].into()),
    }
}

fn render_bond_ref(id: BondId, meta: &impl Metadata) -> Edn<'static> {
    match meta.keyword(Entity::Bond(id)) {
        Some(keyword) => Edn::Keyword(EdnKeyword::owned(keyword.to_string())),
        None => Edn::Int(id.index() as i64),
    }
}

fn render_atom_aliases(meta: &MoleculeMetadata) -> Edn<'static> {
    let aliases = meta.iter_atom_aliases();
    let mut pairs: Vec<Edn<'static>> = Vec::with_capacity(aliases.len() * 2);
    for (name, dsl) in aliases {
        pairs.push(Edn::Keyword(EdnKeyword::owned(name.to_string())));
        pairs.push(dsl.to_edn());
    }
    Edn::Vector(pairs.into())
}

// Unresolved, owned-by-value tree that mirrors the EDN shape. Atom entries and
// per-bond endpoints carry `AtomRef` (index or keyword); constraint leaves carry
// typed per-entity `Constraint*` variants already parsed from their EDN form.
// Lowered destructively via `into_ast(self, cfg)` so that allocations move
// into the final `MoleculeAst`.

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct MoleculeInput {
    pub(crate) atoms: Vec<AtomEntryInput>,
    pub(crate) bonds: Vec<BondEntryInput>,
    pub(crate) dative_bonds: Vec<DativeBondEntryInput>,
    pub(crate) aromatic_systems: Vec<AromaticSystemEntryInput>,
    pub(crate) multicenter_bonds: Vec<MulticenterBondEntryInput>,
    pub(crate) noncovalent_bonds: Vec<NoncovalentBondEntryInput>,
    pub(crate) stereo_atoms: Vec<StereoAtomEntryInput>,
    pub(crate) stereo_bonds: Vec<StereoBondEntryInput>,
    pub(crate) atom_aliases: Vec<(String, Box<AtomDsl>)>,
    pub(crate) constraints: Vec<ConstraintDsl>,
}

/// Atom entry in a parsed molecule map. Mirrors the DSL spec §4 grammar
/// `atom-entry ::= atom-spec | [ keyword atom-spec ]`.
/// TODO: Fix pub(crate) visibility markers on the struct fields.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AtomEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) spec: AtomSpecInput,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum AtomSpecInput {
    Bare(Box<AtomDsl>),
    Alias(String),
}

/// Resolve an atom spec to its `AtomAst`: a bare value is its own atom; an alias is looked up in the
/// table (unknown → error). Shared by the molecule, reaction, and span `into_ast` paths.
pub(super) fn resolve_atom_spec(
    spec: AtomSpecInput,
    namespace: &impl Namespace,
) -> Result<AtomAst, ParseError> {
    match spec {
        AtomSpecInput::Bare(dsl) => Ok(dsl.0),
        AtomSpecInput::Alias(name) => match namespace.find_atom_alias(&name) {
            Some(dsl) => Ok(dsl.0.clone()),
            None => Err(ParseError::InvalidValue(format!(
                "unknown atom alias :{name}"
            ))),
        },
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct BondEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) first: AtomRef,
    pub(crate) second: AtomRef,
    pub(crate) bond: BondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DativeBondEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) donors: Vec<AtomRef>,
    pub(crate) acceptor: AtomRef,
    pub(crate) bond: DativeBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AromaticSystemEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) system: AromaticSystemDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MulticenterBondEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) atoms: Vec<AtomRef>,
    pub(crate) bond: MulticenterBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct NoncovalentBondEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) first: AtomRef,
    pub(crate) second: AtomRef,
    pub(crate) bond: NoncovalentBondDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct StereoAtomEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) site: AtomRef,
    pub(crate) ligands: Vec<StereoLigandRef>,
    pub(crate) stereo: StereoAtomDsl,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct StereoBondEntryInput {
    pub(crate) keyword: Option<String>,
    pub(crate) site: BondRef,
    pub(crate) ligands: Vec<StereoLigandRef>,
    pub(crate) stereo: StereoBondDsl,
}

impl MoleculeInput {
    /// Destructive lowering: consumes the input, resolves refs against the
    /// built id scopes, and produces the final `MoleculeAst` with its
    /// parse-time context. Called from `FromEdn::from_edn` and the streaming path.
    pub(crate) fn into_ast(self) -> Result<(MoleculeAst, MoleculeContext), ParseError> {
        let MoleculeInput {
            atoms: atom_entries,
            bonds: bond_entries,
            dative_bonds: dative_entries,
            aromatic_systems: aromatic_entries,
            multicenter_bonds: multicenter_entries,
            noncovalent_bonds: noncovalent_entries,
            stereo_atoms: stereo_atom_entries,
            stereo_bonds: stereo_bond_entries,
            atom_aliases: alias_entries,
            constraints: constraint_dsls,
        } = self;

        // Register atoms (positions + keywords), then the bijective aliases. `register_*` enforces
        // keyword disjointness (across every entity kind + aliases) and alias bijectivity as it goes, so
        // atom specs resolve against the context once it is complete.
        let mut context = MoleculeContext::default();
        for entry in &atom_entries {
            context.register_atom(entry.keyword.clone())?;
        }
        for (name, dsl) in alias_entries {
            context.register_atom_alias(name, *dsl)?;
        }

        // Resolve atom aliases.
        let atoms: Vec<AtomAst> = atom_entries
            .into_iter()
            .map(|entry| resolve_atom_spec(entry.spec, &context))
            .collect::<Result<_, _>>()?;

        // Bonds.
        let mut bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::with_capacity(bond_entries.len());
        for entry in bond_entries {
            let a = entry.first.resolve(&context)?;
            let b = entry.second.resolve(&context)?;
            context.register_bond(entry.keyword, a, b)?;
            bonds.push((a, b, entry.bond.0));
        }

        // Dative bonds.
        let mut dative_list: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> =
            Vec::with_capacity(dative_entries.len());
        for entry in dative_entries {
            let donors = entry
                .donors
                .into_iter()
                .map(|d| d.resolve(&context))
                .collect::<Result<Vec<_>, _>>()?;
            if donors.is_empty() {
                return Err(ParseError::InvalidValue(
                    "dative bond requires at least one donor".to_string(),
                ));
            }
            let acceptor = entry.acceptor.resolve(&context)?;
            context.register_dative_bond(entry.keyword, &donors, acceptor)?;
            dative_list.push((donors, acceptor, entry.bond.0));
        }

        // Aromatic systems.
        let mut aromatic_list: Vec<(Vec<AtomId>, AromaticSystemAst)> =
            Vec::with_capacity(aromatic_entries.len());
        for entry in aromatic_entries {
            let atoms_resolved: Vec<AtomId> = entry
                .atoms
                .into_iter()
                .map(|r| r.resolve(&context))
                .collect::<Result<_, _>>()?;
            context.register_aromatic_system(entry.keyword, &atoms_resolved)?;
            aromatic_list.push((atoms_resolved, entry.system.0));
        }

        // Multicenter bonds.
        let mut multicenter_list: Vec<(Vec<AtomId>, MulticenterBondAst)> =
            Vec::with_capacity(multicenter_entries.len());
        for entry in multicenter_entries {
            let atoms_resolved: Vec<AtomId> = entry
                .atoms
                .into_iter()
                .map(|r| r.resolve(&context))
                .collect::<Result<_, _>>()?;
            context.register_multicenter_bond(entry.keyword, &atoms_resolved)?;
            multicenter_list.push((atoms_resolved, entry.bond.0));
        }

        // Noncovalent bonds.
        let mut noncovalent_list: Vec<(AtomId, AtomId, NoncovalentBondAst)> =
            Vec::with_capacity(noncovalent_entries.len());
        for entry in noncovalent_entries {
            let first = entry.first.resolve(&context)?;
            let second = entry.second.resolve(&context)?;
            context.register_noncovalent_bond(entry.keyword, first, second)?;
            noncovalent_list.push((first, second, entry.bond.0));
        }

        // Stereo atoms.
        let mut stereo_atom_list: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)> =
            Vec::with_capacity(stereo_atom_entries.len());
        for entry in stereo_atom_entries {
            let site = entry.site.resolve(&context)?;
            let ligands: Vec<StereoLigand> = entry
                .ligands
                .into_iter()
                .map(|l| Ok(StereoLigand::new(l.atom.resolve(&context)?, l.kind)))
                .collect::<Result<_, ParseError>>()?;
            context.register_stereo_atom(entry.keyword, site, &ligands)?;
            stereo_atom_list.push((site, ligands, entry.stereo.0));
        }

        // Stereo bonds.
        let mut stereo_bond_list: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)> =
            Vec::with_capacity(stereo_bond_entries.len());
        for entry in stereo_bond_entries {
            let site = entry.site.resolve(&context)?;
            let ligands: Vec<StereoLigand> = entry
                .ligands
                .into_iter()
                .map(|l| Ok(StereoLigand::new(l.atom.resolve(&context)?, l.kind)))
                .collect::<Result<_, ParseError>>()?;
            context.register_stereo_bond(entry.keyword, site, &ligands)?;
            stereo_bond_list.push((site, ligands, entry.stereo.0));
        }

        // The context is complete; constraints resolve against it directly.
        let constraints = ConstraintsDsl(constraint_dsls).into_ast(&context)?;

        let ast = MoleculeAst::from_parts(MoleculeParts {
            atoms,
            bonds,
            dative: dative_list,
            aromatic: aromatic_list,
            multicenter: multicenter_list,
            noncovalent: noncovalent_list,
            stereo_atoms: stereo_atom_list,
            stereo_bonds: stereo_bond_list,
            constraints,
        });
        Ok((ast, context))
    }
}

pub(super) fn parse_molecule_input(edn: &Edn<'_>) -> Result<MoleculeInput, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "molecule map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let mut input = MoleculeInput::default();
    for (k, v) in m.iter() {
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["molecule".into()],
            });
        };
        match key.name() {
            "atoms" => input.atoms = parse_vec(v, ":atoms", parse_atom_entry)?,
            "bonds" => input.bonds = parse_vec(v, ":bonds", parse_bond_entry)?,
            "dative-bonds" => {
                input.dative_bonds = parse_vec(v, ":dative-bonds", parse_dative_bond_entry)?
            }
            "aromatic-systems" => {
                input.aromatic_systems =
                    parse_vec(v, ":aromatic-systems", parse_aromatic_system_entry)?
            }
            "multicenter-bonds" => {
                input.multicenter_bonds =
                    parse_vec(v, ":multicenter-bonds", parse_multicenter_bond_entry)?
            }
            "noncovalent-bonds" => {
                input.noncovalent_bonds =
                    parse_vec(v, ":noncovalent-bonds", parse_noncovalent_bond_entry)?
            }
            "stereo-atoms" => {
                input.stereo_atoms = parse_vec(v, ":stereo-atoms", parse_stereo_atom_entry)?
            }
            "stereo-bonds" => {
                input.stereo_bonds = parse_vec(v, ":stereo-bonds", parse_stereo_bond_entry)?
            }
            "atom-aliases" => input.atom_aliases = parse_atom_aliases(v)?,
            "constraints" => {
                input.constraints = parse_vec(v, ":constraints", |e| ConstraintDsl::from_edn(e))?
            }
            "guards" => {
                // Spec §4 lists :guards as a future-reserved key; ignore for now.
            }
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["molecule".into()],
                });
            }
        }
    }
    Ok(input)
}

pub(super) fn parse_atom_entry(edn: &Edn<'_>) -> Result<AtomEntryInput, DeError> {
    match edn {
        Edn::Str(s) => {
            let dsl: AtomDsl = s.parse().map_err(|e| DeError::subgrammar("atom", e))?;
            Ok(AtomEntryInput {
                keyword: None,
                spec: AtomSpecInput::Bare(Box::new(dsl)),
            })
        }
        Edn::Keyword(k) => Ok(AtomEntryInput {
            keyword: None,
            spec: AtomSpecInput::Alias(k.name().to_string()),
        }),
        Edn::Vector(v) if v.len() == 2 => {
            let Edn::Keyword(keyword) = &v[0] else {
                return Err(DeError::TypeMismatch {
                    expected: "keyword",
                    got: v[0].kind(),
                    path: vec!["atom-entry".into()],
                });
            };
            let spec = parse_atom_spec(&v[1])?;
            Ok(AtomEntryInput {
                keyword: Some(keyword.name().to_string()),
                spec,
            })
        }
        other => Err(DeError::TypeMismatch {
            expected: "atom-string / keyword / [keyword atom-spec]",
            got: other.kind(),
            path: vec!["atom-entry".into()],
        }),
    }
}

fn parse_atom_spec(edn: &Edn<'_>) -> Result<AtomSpecInput, DeError> {
    match edn {
        Edn::Str(s) => {
            let dsl: AtomDsl = s.parse().map_err(|e| DeError::subgrammar("atom", e))?;
            Ok(AtomSpecInput::Bare(Box::new(dsl)))
        }
        Edn::Keyword(k) => Ok(AtomSpecInput::Alias(k.name().to_string())),
        other => Err(DeError::TypeMismatch {
            expected: "atom-string or keyword alias",
            got: other.kind(),
            path: vec!["atom-spec".into()],
        }),
    }
}

pub(super) fn parse_bond_entry(edn: &Edn<'_>) -> Result<BondEntryInput, DeError> {
    match edn {
        Edn::Vector(v) if v.len() == 3 => Ok(BondEntryInput {
            keyword: None,
            first: AtomRef::from_edn(&v[0])?,
            second: AtomRef::from_edn(&v[1])?,
            bond: BondDsl::from_edn(&v[2])?,
        }),
        Edn::Map(m) => {
            let [a, b] = atoms_pair(m, "bond-entry")?;
            Ok(BondEntryInput {
                keyword: optional_id_keyword(m)?,
                first: a,
                second: b,
                bond: BondDsl::from_edn(required_key(m, "type", "bond-entry")?)?,
            })
        }
        other => Err(DeError::TypeMismatch {
            expected: "bond-entry map or 3-vec",
            got: other.kind(),
            path: vec!["bond-entry".into()],
        }),
    }
}

pub(super) fn parse_dative_bond_entry(edn: &Edn<'_>) -> Result<DativeBondEntryInput, DeError> {
    let m = expect_map(edn, "dative-bond-entry")?;
    let donors = parse_vec(
        required_key(m, "donors", "dative-bond-entry")?,
        ":donors",
        |e| AtomRef::from_edn(e),
    )?;
    Ok(DativeBondEntryInput {
        keyword: optional_id_keyword(m)?,
        donors,
        acceptor: AtomRef::from_edn(required_key(m, "acceptor", "dative-bond-entry")?)?,
        bond: DativeBondDsl::from_edn(required_key(m, "type", "dative-bond-entry")?)?,
    })
}

pub(super) fn parse_aromatic_system_entry(
    edn: &Edn<'_>,
) -> Result<AromaticSystemEntryInput, DeError> {
    let m = expect_map(edn, "aromatic-system-entry")?;
    let system = AromaticSystemDsl::from_edn(required_key(m, "type", "aromatic-system-entry")?)?;
    Ok(AromaticSystemEntryInput {
        keyword: optional_id_keyword(m)?,
        atoms: atoms_vec(m, "aromatic-system-entry")?,
        system,
    })
}

pub(super) fn parse_multicenter_bond_entry(
    edn: &Edn<'_>,
) -> Result<MulticenterBondEntryInput, DeError> {
    let m = expect_map(edn, "multicenter-bond-entry")?;
    let bond = MulticenterBondDsl::from_edn(required_key(m, "type", "multicenter-bond-entry")?)?;
    Ok(MulticenterBondEntryInput {
        keyword: optional_id_keyword(m)?,
        atoms: atoms_vec(m, "multicenter-bond-entry")?,
        bond,
    })
}

pub(super) fn parse_noncovalent_bond_entry(
    edn: &Edn<'_>,
) -> Result<NoncovalentBondEntryInput, DeError> {
    let m = expect_map(edn, "noncovalent-bond-entry")?;
    let [a, b] = atoms_pair(m, "noncovalent-bond-entry")?;
    Ok(NoncovalentBondEntryInput {
        keyword: optional_id_keyword(m)?,
        first: a,
        second: b,
        bond: NoncovalentBondDsl::from_edn(required_key(m, "type", "noncovalent-bond-entry")?)?,
    })
}

pub(super) fn parse_stereo_atom_entry(edn: &Edn<'_>) -> Result<StereoAtomEntryInput, DeError> {
    let m = expect_map(edn, "stereo-atom-entry")?;
    Ok(StereoAtomEntryInput {
        keyword: optional_id_keyword(m)?,
        site: AtomRef::from_edn(required_key(m, "site", "stereo-atom-entry")?)?,
        ligands: parse_vec(
            required_key(m, "ligands", "stereo-atom-entry")?,
            ":ligands",
            parse_stereo_ligand,
        )?,
        stereo: StereoAtomDsl::from_edn(required_key(m, "type", "stereo-atom-entry")?)?,
    })
}

pub(super) fn parse_stereo_bond_entry(edn: &Edn<'_>) -> Result<StereoBondEntryInput, DeError> {
    let m = expect_map(edn, "stereo-bond-entry")?;
    Ok(StereoBondEntryInput {
        keyword: optional_id_keyword(m)?,
        site: BondRef::from_edn(required_key(m, "site", "stereo-bond-entry")?)?,
        ligands: parse_vec(
            required_key(m, "ligands", "stereo-bond-entry")?,
            ":ligands",
            parse_stereo_ligand,
        )?,
        stereo: StereoBondDsl::from_edn(required_key(m, "type", "stereo-bond-entry")?)?,
    })
}

pub(super) fn parse_atom_aliases(edn: &Edn<'_>) -> Result<Vec<(String, Box<AtomDsl>)>, DeError> {
    let Edn::Vector(v) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of keyword/atom-string pairs",
            got: edn.kind(),
            path: vec![":atom-aliases".into()],
        });
    };
    if !v.len().is_multiple_of(2) {
        return Err(DeError::Custom(
            ":atom-aliases must have even length (keyword/atom-string pairs)".into(),
        ));
    }
    let mut out = Vec::with_capacity(v.len() / 2);
    for pair in v.chunks(2) {
        let Edn::Keyword(name) = &pair[0] else {
            return Err(DeError::TypeMismatch {
                expected: "keyword (alias name)",
                got: pair[0].kind(),
                path: vec![":atom-aliases".into()],
            });
        };
        let Edn::Str(s) = &pair[1] else {
            return Err(DeError::TypeMismatch {
                expected: "atom-string",
                got: pair[1].kind(),
                path: vec![":atom-aliases".into()],
            });
        };
        let dsl: AtomDsl = s.parse().map_err(|e| DeError::subgrammar("atom", e))?;
        out.push((name.name().to_string(), Box::new(dsl)));
    }
    Ok(out)
}

fn expect_map<'e>(edn: &'e Edn<'e>, context: &'static str) -> Result<&'e EdnMap<'e>, DeError> {
    match edn {
        Edn::Map(m) => Ok(m),
        other => Err(DeError::TypeMismatch {
            expected: "map",
            got: other.kind(),
            path: vec![context.into()],
        }),
    }
}

#[cfg(test)]
mod tests;
