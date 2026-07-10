//! The `mol!` function-like macro: a compile-checked visual literal that desugars to an L2
//! `MoleculeSpec` and builds a `MoleculeAst`. Grammar (this first slice): comma-separated *paths* of
//! parenthesized atoms joined by bond ops — `(name: elem) - (name: elem) = (other)`, where `elem` is
//! a bare element ident (`C`) or a DSL-spec string (`"C#h3"`); a first mention `(name: elem)` declares,
//! a bare `(name)` references. Undeclared references and duplicate declarations are `compile_error!`s.

use std::collections::HashSet;
use std::iter;

use proc_macro2::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::token::Bracket;
use syn::{bracketed, parenthesized, parse2, Error, Ident, LitStr, Result, Token};

/// The whole `mol!` body: comma-separated paths.
struct MolInput {
    paths: Vec<Path>,
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

/// A bonded chain of atom slots: `first (op atom)*`.
struct Path {
    first: Atom,
    rest: Vec<(Bond, Atom)>,
}

impl Path {
    fn atoms(&self) -> impl Iterator<Item = &Atom> {
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

/// `(name)` (reference) or `(name: element)` (declaration).
struct Atom {
    name: Ident,
    element: Option<ElementSpec>,
}

impl Parse for Atom {
    fn parse(input: ParseStream) -> Result<Self> {
        let content;
        parenthesized!(content in input);
        let name = content.parse::<Ident>()?;
        let element = if content.peek(Token![:]) {
            content.parse::<Token![:]>()?;
            Some(content.parse::<ElementSpec>()?)
        } else {
            None
        };
        Ok(Atom { name, element })
    }
}

/// A bare element ident (`C` → `"C"`) or a DSL-spec string (`"C#h3"`).
enum ElementSpec {
    Bare(Ident),
    Spec(LitStr),
}

impl ElementSpec {
    /// The spec as a string literal for the L2 `Into<AtomAst>` path.
    fn as_lit(&self) -> LitStr {
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

enum Bond {
    Single,
    Double,
    Triple,
    /// `-[ "spec" ]-` — a full DSL bond spec (order, `#a`, charge, spin, ring).
    Spec(LitStr),
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
                let spec = content.parse::<LitStr>()?;
                input.parse::<Token![-]>()?;
                Ok(Bond::Spec(spec))
            } else {
                Ok(Bond::Single)
            }
        }
    }
}

pub fn expand(input: TokenStream) -> TokenStream {
    match parse2::<MolInput>(input).and_then(codegen) {
        Ok(tokens) => tokens,
        Err(error) => error.to_compile_error(),
    }
}

fn codegen(input: MolInput) -> Result<TokenStream> {
    // Pass 1: collect declarations (first mention with `: elem`), reject duplicates.
    let mut declared: Vec<(LitStr, LitStr)> = Vec::new();
    let mut names: HashSet<String> = HashSet::new();
    for path in &input.paths {
        for slot in path.atoms() {
            if let Some(element) = &slot.element {
                if !names.insert(slot.name.to_string()) {
                    return Err(Error::new(
                        slot.name.span(),
                        format!("atom `{}` is declared more than once", slot.name),
                    ));
                }
                let name = LitStr::new(&slot.name.to_string(), slot.name.span());
                declared.push((name, element.as_lit()));
            }
        }
    }
    // Pass 2: every bare reference must resolve to a declaration.
    for path in &input.paths {
        for slot in path.atoms() {
            if slot.element.is_none() && !names.contains(&slot.name.to_string()) {
                return Err(syn::Error::new(
                    slot.name.span(),
                    format!("atom `{}` is referenced but never declared", slot.name),
                ));
            }
        }
    }

    let atoms = if declared.is_empty() {
        quote! {}
    } else {
        let entries = declared
            .iter()
            .map(|(name, element)| quote! { (#name, #element) });
        quote! { + atoms([ #(#entries),* ]) }
    };

    let mut bonds = Vec::new();
    for path in &input.paths {
        let mut previous = &path.first;
        for (op, atom) in &path.rest {
            let a = LitStr::new(&previous.name.to_string(), previous.name.span());
            let b = LitStr::new(&atom.name.to_string(), atom.name.span());
            bonds.push(match op {
                Bond::Single => quote! { single(name(#a), name(#b)) },
                Bond::Double => quote! { double(name(#a), name(#b)) },
                Bond::Triple => quote! { triple(name(#a), name(#b)) },
                Bond::Spec(spec) => quote! { bond(name(#a), name(#b), #spec) },
            });
            previous = atom;
        }
    }

    Ok(quote! {
        {
            #[allow(unused_imports)]
            use ::umol_ast::ast::spec::{atoms, bond, double, name, single, triple, MoleculeSpec};
            (MoleculeSpec::new() #atoms #(+ #bonds)*).build()
        }
    })
}
