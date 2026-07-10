//! Shared `syn` grammar for the `mol!` and `frag!` visual-literal macros: comma-separated *paths* of
//! atoms joined by bond ops, plus the creation-position resolution both macros lower onto. An atom is a
//! named declaration `(name: elem)`, a `(name)` reference, a bare anonymous atom `elem`, or a `^name`
//! port marker (`frag!` only). `elem` is an element ident (`C`) or a DSL-spec string (`"C#h3"`); bond
//! ops are `-` / `=` / `#` / `-[ "spec" ]-`. (`*` is left free for a future Kleene-star operator.)

use std::collections::HashMap;
use std::iter;

use proc_macro2::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::token::{Bracket, Paren};
use syn::{bracketed, parenthesized, Error, Ident, LitStr, Result, Token};

/// The whole macro body: comma-separated paths.
pub(crate) struct MolInput {
    pub(crate) paths: Vec<Path>,
}

impl Parse for MolInput {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut paths = Vec::new();
        while !input.is_empty() {
            paths.push(input.parse::<Path>()?);
            if input.is_empty() {
                break;
            }
            input.parse::<Token![,]>()?;
        }
        Ok(MolInput { paths })
    }
}

/// A bonded chain of atoms: `first (op atom)*`.
pub(crate) struct Path {
    pub(crate) first: Atom,
    pub(crate) rest: Vec<(Bond, Atom)>,
}

impl Path {
    pub(crate) fn atoms(&self) -> impl Iterator<Item = &Atom> {
        iter::once(&self.first).chain(self.rest.iter().map(|(_, atom)| atom))
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
    /// The spec as a string literal for the L2 `Into<AtomAst>` path.
    pub(crate) fn as_lit(&self) -> LitStr {
        match self {
            ElementSpec::Bare(ident) => LitStr::new(&ident.to_string(), ident.span()),
            ElementSpec::Spec(lit) => lit.clone(),
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
    Spec { name: Option<Ident>, spec: LitStr },
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

/// The output of [`resolve_positions`]: the ordered creation specs, and per path each atom
/// occurrence's creation position (`None` marks a `^name` port, which creates no atom).
pub(crate) type Positions = (Vec<LitStr>, Vec<Vec<Option<u32>>>);

/// A label in the shared atom/bond namespace: an atom (with its creation position, for reference
/// resolution) or a bond. Bond labels are inert here — resolving one to its `BondId` for a
/// stereo-overlay reference rides with the overlay work that consumes it.
enum Label {
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
pub(crate) fn resolve_positions(paths: &[Path]) -> Result<Positions> {
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
            if let Bond::Spec { name: Some(name), .. } = op {
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
    Ok((specs, path_positions))
}

/// The L2 bond term for a real-atom-to-real-atom bond, wired by creation position.
pub(crate) fn bond_term(op: &Bond, first: u32, second: u32) -> TokenStream {
    match op {
        Bond::Single => quote! { single(#first, #second) },
        Bond::Double => quote! { double(#first, #second) },
        Bond::Triple => quote! { triple(#first, #second) },
        Bond::Spec { spec, .. } => quote! { bond(#first, #second, #spec) },
    }
}
