//! Shared `syn` grammar for the `mol!` and `frag!` visual-literal macros: comma-separated *paths* of
//! atoms joined by bond ops, plus the creation-position resolution both macros lower onto. An atom is a
//! named declaration `(name: elem)`, a `(name)` reference, a bare anonymous atom `elem`, or a `^name`
//! port marker (`frag!` only). `elem` is an element ident (`C`) or a DSL-spec string (`"C#h3"`); bond
//! ops are `-` / `=` / `#` / `-[ "spec" ]-`. (`*` is left free for a future Kleene-star operator.)

use std::collections::HashMap;

use proc_macro2::{Span, TokenStream};
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::token::{Bracket, Paren};
use syn::{bracketed, parenthesized, Error, Ident, LitStr, Result, Token};

/// Overlay-statement leading keywords.
mod kw {
    syn::custom_keyword!(aromatic);
    syn::custom_keyword!(dative);
    syn::custom_keyword!(multicenter);
    syn::custom_keyword!(noncovalent);
    syn::custom_keyword!(stereo);
    syn::custom_keyword!(atom);
    syn::custom_keyword!(bond);
}

/// The whole macro body: comma-separated statements — a bonded path or an overlay.
pub(crate) struct MolInput {
    pub(crate) statements: Vec<Statement>,
}

impl Parse for MolInput {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut statements = Vec::new();
        while !input.is_empty() {
            statements.push(input.parse::<Statement>()?);
            if input.is_empty() {
                break;
            }
            input.parse::<Token![,]>()?;
        }
        Ok(MolInput { statements })
    }
}

impl MolInput {
    pub(crate) fn paths(&self) -> impl Iterator<Item = &Path> {
        self.statements
            .iter()
            .filter_map(|statement| match statement {
                Statement::Path(path) => Some(path),
                Statement::Overlay(_) => None,
            })
    }

    pub(crate) fn overlays(&self) -> impl Iterator<Item = &Overlay> {
        self.statements
            .iter()
            .filter_map(|statement| match statement {
                Statement::Overlay(overlay) => Some(overlay),
                Statement::Path(_) => None,
            })
    }
}

/// A statement: a bonded path, or an overlay (a keyword-led relation over already-declared atoms/bonds).
pub(crate) enum Statement {
    Path(Path),
    Overlay(Overlay),
}

impl Parse for Statement {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(kw::aromatic)
            || input.peek(kw::dative)
            || input.peek(kw::multicenter)
            || input.peek(kw::noncovalent)
            || input.peek(kw::stereo)
        {
            Ok(Statement::Overlay(input.parse()?))
        } else {
            Ok(Statement::Path(input.parse()?))
        }
    }
}

/// A bonded chain of atoms: `first (op atom)*`.
pub(crate) struct Path {
    pub(crate) first: Atom,
    pub(crate) rest: Vec<(Bond, Atom)>,
}

impl Path {
    pub(crate) fn atoms(&self) -> impl ExactSizeIterator<Item = &Atom> {
        (0..self.rest.len() + 1).map(|index| {
            if index == 0 {
                &self.first
            } else {
                &self.rest[index - 1].1
            }
        })
    }
}

impl Parse for Path {
    fn parse(input: ParseStream) -> Result<Self> {
        let first = input.parse::<Atom>()?;
        let mut rest = Vec::new();
        while input.peek(Token![-]) || input.peek(Token![=]) || input.peek(Token![#]) {
            let op = input.parse::<Bond>()?;
            let atom = input.parse::<Atom>()?;
            rest.push((op, atom));
        }
        Ok(Path { first, rest })
    }
}

/// One atom in a path: a named declaration `(name: elem)`, a `(name)` reference to a declaration, a
/// bare anonymous atom `elem` that nothing can reference, or a `^name` port marker (`frag!` only).
pub(crate) enum Atom {
    Declaration { name: Ident, spec: ElementSpec },
    Reference { name: Ident },
    Anonymous { spec: ElementSpec },
    Port { name: Ident },
}

impl Parse for Atom {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(Token![^]) {
            input.parse::<Token![^]>()?;
            Ok(Atom::Port {
                name: input.parse::<Ident>()?,
            })
        } else if input.peek(Paren) {
            let content;
            parenthesized!(content in input);
            let name = content.parse::<Ident>()?;
            if content.peek(Token![:]) {
                content.parse::<Token![:]>()?;
                Ok(Atom::Declaration {
                    name,
                    spec: content.parse::<ElementSpec>()?,
                })
            } else {
                Ok(Atom::Reference { name })
            }
        } else {
            Ok(Atom::Anonymous {
                spec: input.parse::<ElementSpec>()?,
            })
        }
    }
}

/// A bare element ident (`C` → `"C"`) or a DSL-spec string (`"C#h3"`).
pub(crate) enum ElementSpec {
    Bare(Ident),
    Spec(LitStr),
}

impl ElementSpec {
    /// The spec as a string literal for the L2 `Into<AtomForm>` path.
    pub(crate) fn as_lit(&self) -> LitStr {
        match self {
            ElementSpec::Bare(ident) => LitStr::new(&ident.to_string(), ident.span()),
            ElementSpec::Spec(lit) => lit.clone(),
        }
    }

    pub(crate) fn span(&self) -> Span {
        match self {
            ElementSpec::Bare(ident) => ident.span(),
            ElementSpec::Spec(lit) => lit.span(),
        }
    }
}

impl Parse for ElementSpec {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(LitStr) {
            Ok(ElementSpec::Spec(input.parse()?))
        } else {
            Ok(ElementSpec::Bare(input.parse()?))
        }
    }
}

pub(crate) enum Bond {
    Single,
    Double,
    Triple,
    /// `-[name: "spec"]-` — a full DSL bond spec (order, `#a`, charge, spin, ring), with an optional
    /// label binding the bond in the shared atom/bond namespace.
    Spec {
        name: Option<Ident>,
        spec: LitStr,
    },
}

impl Parse for Bond {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(Token![=]) {
            input.parse::<Token![=]>()?;
            Ok(Bond::Double)
        } else if input.peek(Token![#]) {
            input.parse::<Token![#]>()?;
            Ok(Bond::Triple)
        } else {
            input.parse::<Token![-]>()?;
            if input.peek(Bracket) {
                let content;
                bracketed!(content in input);
                let name = if content.peek(Ident) && content.peek2(Token![:]) {
                    let name = content.parse::<Ident>()?;
                    content.parse::<Token![:]>()?;
                    Some(name)
                } else {
                    None
                };
                let spec = content.parse::<LitStr>()?;
                input.parse::<Token![-]>()?;
                Ok(Bond::Spec { name, spec })
            } else {
                Ok(Bond::Single)
            }
        }
    }
}

/// A ligand in a stereo overlay: an atom reference, or a virtual `"#h"` / `"#n"` placeholder.
pub(crate) enum Ligand {
    Atom(Atom),
    ImplicitHydrogen,
    LonePair,
}

impl Parse for Ligand {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(LitStr) {
            let placeholder: LitStr = input.parse()?;
            match placeholder.value().as_str() {
                "#h" => Ok(Ligand::ImplicitHydrogen),
                "#n" => Ok(Ligand::LonePair),
                other => Err(Error::new(
                    placeholder.span(),
                    format!("unknown ligand placeholder {other:?}; use \"#h\" or \"#n\""),
                )),
            }
        } else {
            Ok(Ligand::Atom(input.parse()?))
        }
    }
}

/// An overlay statement: a keyword-led relation over already-declared atoms/bonds (references only).
/// Each desugars to its L2 term; the optional payload is a quoted DSL spec for the entity's own Ast.
pub(crate) enum Overlay {
    Aromatic {
        atoms: Vec<Atom>,
        payload: Option<LitStr>,
    },
    Dative {
        donors: Vec<Atom>,
        acceptor: Atom,
        payload: Option<LitStr>,
    },
    Multicenter {
        atoms: Vec<Atom>,
        payload: Option<LitStr>,
    },
    Noncovalent {
        atoms: Vec<Atom>,
        payload: Option<LitStr>,
    },
    StereoAtom {
        site: Atom,
        ligands: Vec<Ligand>,
        payload: Option<LitStr>,
    },
    StereoBond {
        site: Atom,
        ligands: Vec<Ligand>,
        payload: Option<LitStr>,
    },
}

impl Parse for Overlay {
    fn parse(input: ParseStream) -> Result<Self> {
        if input.peek(kw::aromatic) {
            input.parse::<kw::aromatic>()?;
            let atoms = bracketed_list(input)?;
            Ok(Overlay::Aromatic {
                atoms,
                payload: payload(input)?,
            })
        } else if input.peek(kw::dative) {
            input.parse::<kw::dative>()?;
            let donors = bracketed_list(input)?;
            let acceptor = input.parse::<Atom>()?;
            Ok(Overlay::Dative {
                donors,
                acceptor,
                payload: payload(input)?,
            })
        } else if input.peek(kw::multicenter) {
            input.parse::<kw::multicenter>()?;
            let atoms = bracketed_list(input)?;
            Ok(Overlay::Multicenter {
                atoms,
                payload: payload(input)?,
            })
        } else if input.peek(kw::noncovalent) {
            input.parse::<kw::noncovalent>()?;
            let atoms = bracketed_list(input)?;
            Ok(Overlay::Noncovalent {
                atoms,
                payload: payload(input)?,
            })
        } else {
            input.parse::<kw::stereo>()?;
            if input.peek(kw::atom) {
                input.parse::<kw::atom>()?;
                let site = input.parse::<Atom>()?;
                let ligands = bracketed_list(input)?;
                Ok(Overlay::StereoAtom {
                    site,
                    ligands,
                    payload: payload(input)?,
                })
            } else {
                input.parse::<kw::bond>()?;
                let site = input.parse::<Atom>()?;
                let ligands = bracketed_list(input)?;
                Ok(Overlay::StereoBond {
                    site,
                    ligands,
                    payload: payload(input)?,
                })
            }
        }
    }
}

/// Parse a whitespace-separated `[ … ]` list.
fn bracketed_list<T: Parse>(input: ParseStream) -> Result<Vec<T>> {
    let content;
    bracketed!(content in input);
    let mut items = Vec::new();
    while !content.is_empty() {
        items.push(content.parse()?);
    }
    Ok(items)
}

/// Parse an optional `: "payload"` trailer.
fn payload(input: ParseStream) -> Result<Option<LitStr>> {
    if input.peek(Token![:]) {
        input.parse::<Token![:]>()?;
        Ok(Some(input.parse()?))
    } else {
        Ok(None)
    }
}

/// The output of [`resolve_positions`]: the ordered creation specs, per-path atom positions (`None`
/// marks a `^name` port, which creates no atom), and the shared label namespace (for overlay refs).
pub(crate) struct Resolved {
    pub(crate) specs: Vec<LitStr>,
    pub(crate) path_positions: Vec<Vec<Option<u32>>>,
    pub(crate) labels: HashMap<String, Label>,
}

/// A label in the shared atom/bond namespace: an atom (with its creation position, for reference
/// resolution) or a bond (referenced by name — the stereo-overlay site).
pub(crate) enum Label {
    Atom(u32),
    Bond,
}

/// Register `name` in the shared namespace; a duplicate — atom or bond — is an error, since atom and
/// bond labels share one namespace.
fn insert_label(labels: &mut HashMap<String, Label>, name: &Ident, label: Label) -> Result<()> {
    if let Some(existing) = labels.insert(name.to_string(), label) {
        let existing_kind = match existing {
            Label::Atom(_) => "an atom",
            Label::Bond => "a bond",
        };
        return Err(Error::new(
            name.span(),
            format!("label `{name}` is already declared as {existing_kind}; atom and bond labels share one namespace"),
        ));
    }
    Ok(())
}

/// Resolve an atom reference `(name)` to its creation position. A name bound to a bond, or unknown, is
/// an error.
fn resolve_atom_ref(labels: &HashMap<String, Label>, name: &Ident) -> Result<u32> {
    match labels.get(&name.to_string()) {
        Some(Label::Atom(position)) => Ok(*position),
        Some(Label::Bond) => Err(Error::new(
            name.span(),
            format!("`{name}` is a bond label, not an atom"),
        )),
        None => Err(Error::new(
            name.span(),
            format!("atom `{name}` is referenced but never declared"),
        )),
    }
}

/// Assign a creation position to every declaration and anonymous atom in appearance order, collecting
/// their specs and registering declaration and `-[name: …]-` bond labels in the shared namespace
/// (duplicates rejected). Then resolve each atom occurrence to its position — declarations and
/// anonymous atoms advance the counter in that same order, references resolve to their declaration,
/// and `^name` port markers resolve to `None` (they create no atom).
pub(crate) fn resolve_positions(paths: &[&Path]) -> Result<Resolved> {
    let mut labels: HashMap<String, Label> = HashMap::new();
    let mut specs: Vec<LitStr> = Vec::new();
    for path in paths {
        for atom in path.atoms() {
            match atom {
                Atom::Declaration { name, spec } => {
                    insert_label(&mut labels, name, Label::Atom(specs.len() as u32))?;
                    specs.push(spec.as_lit());
                }
                Atom::Anonymous { spec } => specs.push(spec.as_lit()),
                Atom::Reference { .. } | Atom::Port { .. } => {}
            }
        }
    }
    for path in paths {
        for (op, _) in &path.rest {
            if let Bond::Spec {
                name: Some(name), ..
            } = op
            {
                insert_label(&mut labels, name, Label::Bond)?;
            }
        }
    }

    let mut next_position = 0u32;
    let mut path_positions: Vec<Vec<Option<u32>>> = Vec::new();
    for path in paths {
        let mut row = Vec::new();
        for atom in path.atoms() {
            let position = match atom {
                Atom::Declaration { .. } | Atom::Anonymous { .. } => {
                    let position = next_position;
                    next_position += 1;
                    Some(position)
                }
                Atom::Reference { name } => Some(resolve_atom_ref(&labels, name)?),
                Atom::Port { .. } => None,
            };
            row.push(position);
        }
        path_positions.push(row);
    }
    Ok(Resolved {
        specs,
        path_positions,
        labels,
    })
}

/// Resolve a bond reference `(name)` (a stereo-overlay site). A name bound to an atom, or unknown, is
/// an error.
fn resolve_bond_ref(labels: &HashMap<String, Label>, name: &Ident) -> Result<()> {
    match labels.get(&name.to_string()) {
        Some(Label::Bond) => Ok(()),
        Some(Label::Atom(_)) => Err(Error::new(
            name.span(),
            format!("`{name}` is an atom, not a bond"),
        )),
        None => Err(Error::new(
            name.span(),
            format!("bond `{name}` is referenced but never declared"),
        )),
    }
}

/// The L2 bond term for a real-atom-to-real-atom bond, wired by creation position. A `-[name: …]-`
/// bond emits `named_bond` so a later stereo-bond site can reference it by name.
pub(crate) fn bond_term(op: &Bond, first: u32, second: u32) -> TokenStream {
    match op {
        Bond::Single => quote! { single(#first, #second) },
        Bond::Double => quote! { double(#first, #second) },
        Bond::Triple => quote! { triple(#first, #second) },
        Bond::Spec {
            name: Some(name),
            spec,
        } => {
            let name = name.to_string();
            quote! { named_bond(#name, #first, #second, #spec) }
        }
        Bond::Spec { name: None, spec } => quote! { bond(#first, #second, #spec) },
    }
}

/// An overlay atom reference → an `AtomArg::Index` position. Overlays reference declared atoms only;
/// an inline declaration or anonymous atom in an overlay is an error (declare it in a path).
fn overlay_atom(atom: &Atom, labels: &HashMap<String, Label>) -> Result<TokenStream> {
    match atom {
        Atom::Reference { name } => {
            let position = resolve_atom_ref(labels, name)?;
            Ok(quote! { AtomArg::Index(#position) })
        }
        Atom::Declaration { name, .. } => Err(Error::new(
            name.span(),
            "overlay participants reference declared atoms — declare atoms in a path",
        )),
        Atom::Anonymous { spec } => Err(Error::new(
            spec.span(),
            "overlay participants reference declared atoms — declare atoms in a path",
        )),
        Atom::Port { name } => Err(Error::new(
            name.span(),
            "a port cannot appear in an overlay",
        )),
    }
}

/// The `ast: impl Into<…>` payload argument of an overlay term — the quoted DSL spec, or the entity's
/// `default()` when none was given.
fn overlay_payload(payload: &Option<LitStr>, ast_ty: TokenStream) -> TokenStream {
    match payload {
        Some(spec) => quote! { #spec },
        None => quote! { <#ast_ty>::default() },
    }
}

/// The L2 term for one overlay, resolving its references against the shared label namespace.
pub(crate) fn overlay_term(
    overlay: &Overlay,
    labels: &HashMap<String, Label>,
) -> Result<TokenStream> {
    match overlay {
        Overlay::Aromatic { atoms, payload } => {
            let atoms = overlay_atoms(atoms, labels)?;
            let ast = overlay_payload(payload, quote! { ::umol_graph_ir::ir::AromaticSystemForm });
            Ok(quote! { aromatic_system([ #(#atoms),* ], #ast) })
        }
        Overlay::Dative {
            donors,
            acceptor,
            payload,
        } => {
            let donors = overlay_atoms(donors, labels)?;
            let acceptor = overlay_atom(acceptor, labels)?;
            let ast = overlay_payload(payload, quote! { ::umol_graph_ir::ir::DativeBondForm });
            Ok(quote! { dative_bond([ #(#donors),* ], #acceptor, #ast) })
        }
        Overlay::Multicenter { atoms, payload } => {
            let atoms = overlay_atoms(atoms, labels)?;
            let ast = overlay_payload(payload, quote! { ::umol_graph_ir::ir::MulticenterBondForm });
            Ok(quote! { multicenter_bond([ #(#atoms),* ], #ast) })
        }
        Overlay::Noncovalent { atoms, payload } => {
            if atoms.len() != 2 {
                return Err(Error::new(
                    Span::call_site(),
                    "noncovalent takes exactly two atoms",
                ));
            }
            let first = overlay_atom(&atoms[0], labels)?;
            let second = overlay_atom(&atoms[1], labels)?;
            let ast = overlay_payload(payload, quote! { ::umol_graph_ir::ir::NoncovalentBondForm });
            Ok(quote! { noncovalent_bond(#first, #second, #ast) })
        }
        Overlay::StereoAtom {
            site,
            ligands,
            payload,
        } => {
            let site = overlay_atom(site, labels)?;
            let ligands = overlay_ligands(ligands, labels)?;
            let ast = overlay_payload(payload, quote! { ::umol_graph_ir::ir::StereoAtomAst });
            Ok(quote! { stereo_atom(#site, [ #(#ligands),* ], #ast) })
        }
        Overlay::StereoBond {
            site,
            ligands,
            payload,
        } => {
            let Atom::Reference { name } = site else {
                return Err(Error::new(
                    Span::call_site(),
                    "a stereo-bond site must reference a named bond",
                ));
            };
            resolve_bond_ref(labels, name)?;
            let site = name.to_string();
            let ligands = overlay_ligands(ligands, labels)?;
            let ast = overlay_payload(payload, quote! { ::umol_graph_ir::ir::StereoBondAst });
            Ok(quote! { stereo_bond(#site, [ #(#ligands),* ], #ast) })
        }
    }
}

fn overlay_atoms(atoms: &[Atom], labels: &HashMap<String, Label>) -> Result<Vec<TokenStream>> {
    atoms
        .iter()
        .map(|atom| overlay_atom(atom, labels))
        .collect()
}

fn overlay_ligands(
    ligands: &[Ligand],
    labels: &HashMap<String, Label>,
) -> Result<Vec<TokenStream>> {
    ligands
        .iter()
        .map(|ligand| match ligand {
            Ligand::Atom(atom) => {
                let atom = overlay_atom(atom, labels)?;
                Ok(quote! { StereoLigandArg::Atom(#atom) })
            }
            Ligand::ImplicitHydrogen => Ok(quote! { StereoLigandArg::ImplicitHydrogen }),
            Ligand::LonePair => Ok(quote! { StereoLigandArg::LonePair }),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::single("C", vec!["C"])]
    #[case::chain("C-O=N", vec!["C", "O", "N"])]
    fn test_path_atoms(#[case] input: &str, #[case] expected: Vec<&str>) {
        let path = syn::parse_str::<Path>(input).expect("path fixture must parse");
        let mut atoms = path.atoms();
        assert_eq!(atoms.len(), expected.len());
        assert_eq!(atoms.size_hint(), (expected.len(), Some(expected.len())),);
        for expected_atom in expected {
            let previous = atoms.len();
            let actual = atoms.next().map(|atom| match atom {
                Atom::Declaration { spec, .. } | Atom::Anonymous { spec } => spec.as_lit().value(),
                Atom::Reference { name } | Atom::Port { name } => name.to_string(),
            });
            assert_eq!(actual.as_deref(), Some(expected_atom));
            let remaining = atoms.len();
            assert_eq!(remaining, previous - 1);
            assert_eq!(atoms.size_hint(), (remaining, Some(remaining)));
        }
        assert_eq!(atoms.next().map(|_| ()), None);
        assert_eq!(atoms.len(), 0);
    }
}
