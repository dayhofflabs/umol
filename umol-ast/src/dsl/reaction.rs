//! Reaction DSL.
//!
//! `ReactionDsl` wraps a `ReactionAst` together with the `ReactionMetadata` that records
//! the surface-form bindings (the lhs molecule metadata plus created-entity id ↔ name and
//! atom-alias bindings). The EDN form is a map keyed by `:lhs` (a molecule map) and
//! `:deltas` (a vector of `:add` / `:remove` / `:modify` / `:constraint` operations). Each
//! entity delegates to its own entity DSL.

use bimap::BiBTreeMap;
use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn};

use super::atom::{lower_atom, raise_atom, AtomDsl, PartialAtomDsl};
use super::bond::{lower_bond, raise_bond, PartialBondDsl};
use super::config::{DeltaDefaults, ReactionDefaults};
use super::constraint::{read_constraint_dsl, ConstraintDsl};
use super::edn_utils::{
    consume_single_key_map_close, parse_single_key_map, read_single_key_map_header,
};
use super::molecule::{
    parse_atom_entry, parse_bond_entry, read_atom_entry, read_bond_entry, AtomEntryInput,
    BondEntryInput, MoleculeDsl, MoleculeInput, MoleculeMetadata,
};
use super::refs::{read_atom_ref, read_bond_ref, AtomRef, BondRef};
use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::delta::{AtomDelta, BondDelta, Delta};
use crate::ast::id::{AtomId, BondId};
use crate::ast::reaction::ReactionAst;
use crate::ast::traits::{FromAst, IntoAst};

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
        _ => {}
    }
}

/// Surface-form metadata paired with a `ReactionAst`: the lhs molecule metadata plus the
/// created-entity id bindings and atom aliases introduced by the deltas. Mirrors
/// `MoleculeMetadata` for the atom/bond entities (the reaction admits the `[:C "C#h3"]`
/// alias notation for added atoms).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReactionMetadata {
    lhs: MoleculeMetadata,
    atom_ids: IndexMap<AtomId, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
    bond_ids: IndexMap<BondId, String>,
}

impl ReactionMetadata {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn lhs(&self) -> &MoleculeMetadata {
        &self.lhs
    }

    pub fn atom_id(&self, id: AtomId) -> Option<&str> {
        self.atom_ids.get(&id).map(String::as_str)
    }

    pub fn bond_id(&self, id: BondId) -> Option<&str> {
        self.bond_ids.get(&id).map(String::as_str)
    }

    /// Name of the alias bound to this atom DSL, if any.
    pub fn atom_alias_for(&self, dsl: &AtomDsl) -> Option<&str> {
        self.atom_aliases.get_by_right(dsl).map(String::as_str)
    }

    pub fn has_atom_alias(&self, name: &str) -> bool {
        self.atom_aliases.contains_left(name)
    }

    pub fn has_atom_aliases(&self) -> bool {
        !self.atom_aliases.is_empty()
    }

    pub fn atom_aliases_len(&self) -> usize {
        self.atom_aliases.len()
    }

    pub fn iter_atom_aliases(&self) -> impl Iterator<Item = (&str, &AtomDsl)> {
        self.atom_aliases
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_ref()))
    }

    pub fn set_atom_id(&mut self, id: AtomId, name: impl Into<String>) {
        self.atom_ids.insert(id, name.into());
    }

    pub fn set_bond_id(&mut self, id: BondId, name: impl Into<String>) {
        self.bond_ids.insert(id, name.into());
    }

    /// Insert an atom alias. Last-wins on either side of the bijection: a
    /// duplicate name displaces its prior atom-dsl mapping, and a duplicate
    /// atom-dsl displaces its prior name. Callers that need collision
    /// detection check upstream.
    pub fn add_atom_alias(&mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) {
        self.atom_aliases.insert(name.into(), Box::new(atom.into()));
    }

    pub fn with_atom_id(mut self, id: AtomId, name: impl Into<String>) -> Self {
        self.set_atom_id(id, name);
        self
    }

    pub fn with_bond_id(mut self, id: BondId, name: impl Into<String>) -> Self {
        self.set_bond_id(id, name);
        self
    }

    pub fn with_atom_alias(mut self, name: impl Into<String>, atom: impl Into<AtomDsl>) -> Self {
        self.add_atom_alias(name, atom);
        self
    }
}

/// One unresolved delta parsed from a `:deltas` entry. Refs stay symbolic
/// (`AtomRef`/`BondRef`) and an `:add` carries the molecule entry verbatim
/// (bare atom / alias / created-id resolved in R7); a `:modify` RHS is a
/// partial entity AST (unspecified fields `Undetermined`).
#[derive(Debug, PartialEq)]
pub(crate) enum DeltaInput {
    AtomAdd(AtomEntryInput),
    AtomRemove(AtomRef),
    AtomModify(AtomRef, AtomAst),
    BondAdd(BondEntryInput),
    BondRemove(BondRef),
    BondModify(BondRef, BondAst),
    ConstraintAdd(ConstraintDsl),
    ConstraintRemove(ConstraintDsl),
}

/// Raw parse target for a reaction: the lhs molecule input plus the unresolved
/// deltas. Resolution (`into_ast`, R7) lifts this to `(ReactionAst, ReactionMetadata)`.
#[derive(Debug)]
pub(crate) struct ReactionInput {
    lhs: MoleculeInput,
    deltas: Vec<DeltaInput>,
}

fn read_delta_input(de: &mut EdnStreamDeserializer<'_>) -> Result<DeltaInput, EdnError> {
    let entity = read_single_key_map_header(de)?;
    let input = match entity.as_str() {
        "atom" => read_delta_atom_input(de)?,
        "bond" => read_delta_bond_input(de)?,
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
            let dsl: PartialAtomDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("partial-atom", e))?;
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
            let dsl: PartialBondDsl = s
                .as_ref()
                .parse()
                .map_err(|e| DeError::subgrammar("partial-bond", e))?;
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
                PartialAtomDsl::from_edn(&v[1])?.0,
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
                PartialBondDsl::from_edn(&v[1])?.0,
            ))
        }
        o => Err(DeError::Custom(format!("unknown bond delta op :{o}"))),
    }
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::atom::ElementAst;
    use crate::ast::constraint::Constraint;
    use crate::ast::delta::{ConstraintDelta, Deltas};
    use crate::ast::edit::AtomFieldChange;
    use crate::ast::value::ValueAst;
    use crate::dsl::bond::BondDsl;
    use crate::dsl::constraint::MoleculeConstraintDsl;
    use crate::dsl::molecule::AtomSpecInput;
    use crate::mol;

    #[rstest]
    #[case::sn2(ReactionAst::new(
        mol!(r##"{:atoms ["C" "Br"] :bonds [[0 1 "1"]]}"##),
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
            id: None,
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                a.implicit_hydrogens = ValueAst::Lit(3);
                a
            }))),
        })
    )]
    #[case::add_id(
        r##"{:atom {:add [:nu "O#h1"]}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            id: Some("nu".into()),
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
            id: None,
            spec: AtomSpecInput::Alias("foo".into()),
        })
    )]
    #[case::remove_id("{:atom {:remove :br}}", DeltaInput::AtomRemove(AtomRef::Id("br".into())))]
    #[case::remove_index("{:atom {:remove 1}}", DeltaInput::AtomRemove(AtomRef::Index(1)))]
    #[case::modify(
        r##"{:atom {:modify [:br "#c-1"]}}"##,
        DeltaInput::AtomModify(AtomRef::Id("br".into()), {
            let mut a = AtomAst::new(ElementAst::Undetermined);
            a.charge = ValueAst::Lit(-1);
            a
        })
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
            id: None,
            spec: AtomSpecInput::Bare(Box::new(AtomDsl({
                let mut a = AtomAst::new(ElementAst::Lit(Element::C));
                a.implicit_hydrogens = ValueAst::Lit(3);
                a
            }))),
        })
    )]
    #[case::add_id(
        r##"{:atom {:add [:nu "O#h1"]}}"##,
        DeltaInput::AtomAdd(AtomEntryInput {
            id: Some("nu".into()),
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
            id: None,
            spec: AtomSpecInput::Alias("foo".into()),
        })
    )]
    #[case::remove_id("{:atom {:remove :br}}", DeltaInput::AtomRemove(AtomRef::Id("br".into())))]
    #[case::remove_index("{:atom {:remove 1}}", DeltaInput::AtomRemove(AtomRef::Index(1)))]
    #[case::modify(
        r##"{:atom {:modify [:br "#c-1"]}}"##,
        DeltaInput::AtomModify(AtomRef::Id("br".into()), {
            let mut a = AtomAst::new(ElementAst::Undetermined);
            a.charge = ValueAst::Lit(-1);
            a
        })
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
            id: None,
            a: AtomRef::Index(0),
            b: AtomRef::Index(1),
            bond: BondDsl(BondAst::from_order(1)),
        })
    )]
    #[case::add_map_id(
        r##"{:bond {:add {:id :b1 :atoms [:c :nu] :type "2"}}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            id: Some("b1".into()),
            a: AtomRef::Id("c".into()),
            b: AtomRef::Id("nu".into()),
            bond: BondDsl(BondAst::from_order(2)),
        })
    )]
    #[case::remove_id("{:bond {:remove :b1}}", DeltaInput::BondRemove(BondRef::Id("b1".into())))]
    #[case::remove_index("{:bond {:remove 0}}", DeltaInput::BondRemove(BondRef::Index(0)))]
    #[case::modify(
        r##"{:bond {:modify [:b1 "2"]}}"##,
        DeltaInput::BondModify(BondRef::Id("b1".into()), BondAst::from_order(2))
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
            id: None,
            a: AtomRef::Index(0),
            b: AtomRef::Index(1),
            bond: BondDsl(BondAst::from_order(1)),
        })
    )]
    #[case::add_map_id(
        r##"{:bond {:add {:id :b1 :atoms [:c :nu] :type "2"}}}"##,
        DeltaInput::BondAdd(BondEntryInput {
            id: Some("b1".into()),
            a: AtomRef::Id("c".into()),
            b: AtomRef::Id("nu".into()),
            bond: BondDsl(BondAst::from_order(2)),
        })
    )]
    #[case::remove_id("{:bond {:remove :b1}}", DeltaInput::BondRemove(BondRef::Id("b1".into())))]
    #[case::remove_index("{:bond {:remove 0}}", DeltaInput::BondRemove(BondRef::Index(0)))]
    #[case::modify(
        r##"{:bond {:modify [:b1 "2"]}}"##,
        DeltaInput::BondModify(BondRef::Id("b1".into()), BondAst::from_order(2))
    )]
    fn test_read_delta_input_bond(#[case] input: &str, #[case] expected: DeltaInput) {
        assert_eq!(
            read_delta_input(&mut EdnStreamDeserializer::new(input)).unwrap(),
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
}
