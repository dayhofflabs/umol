//! The `mol!` function-like macro: a compile-checked visual literal that desugars to an L2
//! `MoleculeSpec` and builds a `MoleculeAst`. Grammar (this first slice): comma-separated *paths* of
//! parenthesized atoms joined by bond ops — `(name: elem) - (name: elem) = (other)`, where `elem` is
//! a bare element ident (`C`) or a DSL-spec string (`"C#h3"`); a first mention `(name: elem)` declares,
//! a bare `(name)` references. Undeclared references and duplicate declarations are `compile_error!`s.

use std::collections::HashSet;

use proc_macro2::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::{parenthesized, Ident, LitStr, Token};

/// The whole `mol!` body: comma-separated paths.
struct MolInput {
    paths: Vec<Path>,
}

impl Parse for MolInput {
    fn parse(input: ParseStream) -> syn::Result<Self> {
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
    first: AtomSlot,
    rest: Vec<(BondOp, AtomSlot)>,
}

impl Path {
    fn slots(&self) -> impl Iterator<Item = &AtomSlot> {
        std::iter::once(&self.first).chain(self.rest.iter().map(|(_, atom)| atom))
    }
}

impl Parse for Path {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        let first = input.parse::<AtomSlot>()?;
        let mut rest = Vec::new();
        while input.peek(Token![-]) || input.peek(Token![=]) || input.peek(Token![#]) {
            let op = input.parse::<BondOp>()?;
            let atom = input.parse::<AtomSlot>()?;
            rest.push((op, atom));
        }
        Ok(Path { first, rest })
    }
}

/// `(name)` (reference) or `(name: element)` (declaration).
struct AtomSlot {
    name: Ident,
    element: Option<ElementSpec>,
}

impl Parse for AtomSlot {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        let content;
        parenthesized!(content in input);
        let name = content.parse::<Ident>()?;
        let element = if content.peek(Token![:]) {
            content.parse::<Token![:]>()?;
            Some(content.parse::<ElementSpec>()?)
        } else {
            None
        };
        Ok(AtomSlot { name, element })
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
    fn parse(input: ParseStream) -> syn::Result<Self> {
        if input.peek(LitStr) {
            Ok(ElementSpec::Spec(input.parse()?))
        } else {
            Ok(ElementSpec::Bare(input.parse()?))
        }
    }
}

enum BondOp {
    Single,
    Double,
    Triple,
}

impl Parse for BondOp {
    fn parse(input: ParseStream) -> syn::Result<Self> {
        if input.peek(Token![=]) {
            input.parse::<Token![=]>()?;
            Ok(BondOp::Double)
        } else if input.peek(Token![#]) {
            input.parse::<Token![#]>()?;
            Ok(BondOp::Triple)
        } else {
            input.parse::<Token![-]>()?;
            Ok(BondOp::Single)
        }
    }
}

pub fn expand(input: TokenStream) -> TokenStream {
    match syn::parse2::<MolInput>(input).and_then(codegen) {
        Ok(tokens) => tokens,
        Err(error) => error.to_compile_error(),
    }
}

fn codegen(input: MolInput) -> syn::Result<TokenStream> {
    // Pass 1: collect declarations (first mention with `: elem`), reject duplicates.
    let mut declared: Vec<(LitStr, LitStr)> = Vec::new();
    let mut names: HashSet<String> = HashSet::new();
    for path in &input.paths {
        for slot in path.slots() {
            if let Some(element) = &slot.element {
                if !names.insert(slot.name.to_string()) {
                    return Err(syn::Error::new(
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
        for slot in path.slots() {
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
        let entries = declared.iter().map(|(name, element)| quote! { (#name, #element) });
        quote! { + atoms([ #(#entries),* ]) }
    };

    let mut bonds = Vec::new();
    for path in &input.paths {
        let mut previous = &path.first;
        for (op, atom) in &path.rest {
            let a = LitStr::new(&previous.name.to_string(), previous.name.span());
            let b = LitStr::new(&atom.name.to_string(), atom.name.span());
            bonds.push(match op {
                BondOp::Single => quote! { single(name(#a), name(#b)) },
                BondOp::Double => quote! { double(name(#a), name(#b)) },
                BondOp::Triple => quote! { triple(name(#a), name(#b)) },
            });
            previous = atom;
        }
    }

    Ok(quote! {
        {
            #[allow(unused_imports)]
            use ::umol_ast::ast::spec::{atoms, double, name, single, triple, MoleculeSpec};
            (MoleculeSpec::new() #atoms #(+ #bonds)*).build()
        }
    })
}
