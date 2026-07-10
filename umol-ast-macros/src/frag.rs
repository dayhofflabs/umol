//! The `frag!` function-like macro: the shared path grammar (see [`crate::parse`]) extended with
//! `^name` port markers, desugaring to an L3 `Fragment` — a `MoleculeAst` body built via L2, wrapped
//! with one port per marker. A bond incident to a `^name` marker declares a port on the real endpoint
//! (name `name`, colour = the bond's `BondAst`); a bond between two markers is a `compile_error!`.

use proc_macro2::TokenStream;
use quote::quote;
use syn::{parse2, Error, Ident, LitStr, Result};

use crate::parse::{bond_term, resolve_positions, Atom, Bond, MolInput};

pub fn expand(input: TokenStream) -> TokenStream {
    match parse2::<MolInput>(input).and_then(codegen) {
        Ok(tokens) => tokens,
        Err(error) => error.to_compile_error(),
    }
}

fn codegen(input: MolInput) -> Result<TokenStream> {
    let (specs, positions) = resolve_positions(&input.paths)?;

    let atoms = if specs.is_empty() {
        quote! {}
    } else {
        quote! { + atoms([ #(#specs),* ]) }
    };

    let mut bonds = Vec::new();
    let mut port_calls = Vec::new();
    for (path, row) in input.paths.iter().zip(&positions) {
        let occurrences: Vec<&Atom> = path.atoms().collect();
        for (index, (op, _)) in path.rest.iter().enumerate() {
            match (row[index], row[index + 1]) {
                (Some(first), Some(second)) => bonds.push(bond_term(op, first, second)),
                (Some(atom), None) => {
                    port_calls.push(port_call(port_name(occurrences[index + 1]), atom, op));
                }
                (None, Some(atom)) => {
                    port_calls.push(port_call(port_name(occurrences[index]), atom, op));
                }
                (None, None) => {
                    return Err(Error::new(
                        port_name(occurrences[index]).span(),
                        "a bond cannot connect two ports",
                    ));
                }
            }
        }
    }

    Ok(quote! {
        {
            #[allow(unused_imports)]
            use ::umol_ast::ast::spec::{atoms, bond, double, single, triple, MoleculeSpec};
            #[allow(unused_imports)]
            use ::umol_ast::ast::{AtomId, BondAst, Fragment};
            Fragment::new((MoleculeSpec::new() #atoms #(+ #bonds)*).build())
                #(#port_calls)*
        }
    })
}

/// The `^name` marker's ident. A `None` position (per [`resolve_positions`]) is always a port marker.
fn port_name(atom: &Atom) -> &Ident {
    match atom {
        Atom::Port { name } => name,
        _ => unreachable!("a `None` position is always a port marker"),
    }
}

/// A `.with_port(name, AtomId(atom), colour)` call — the port's colour is the bond op's `BondAst`.
fn port_call(name: &Ident, atom: u32, op: &Bond) -> TokenStream {
    let name = LitStr::new(&name.to_string(), name.span());
    let colour = match op {
        Bond::Single => quote! { BondAst::from_order(1) },
        Bond::Double => quote! { BondAst::from_order(2) },
        Bond::Triple => quote! { BondAst::from_order(3) },
        Bond::Spec(spec) => quote! { #spec },
    };
    quote! { .with_port(#name, AtomId(#atom), #colour) }
}
