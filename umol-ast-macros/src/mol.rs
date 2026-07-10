//! The `mol!` function-like macro: a compile-checked visual literal over the shared path grammar
//! (see [`crate::parse`]) that desugars to an L2 `MoleculeSpec` and builds a `MoleculeAst`. Every atom
//! resolves to a creation position; the molecule is emitted as one nameless `atoms([spec, …])` term
//! wired by position. Port markers (`^name`) belong to `frag!` and are rejected here.

use proc_macro2::TokenStream;
use quote::quote;
use syn::{parse2, Error, Result};

use crate::parse::{bond_term, resolve_positions, Atom, MolInput};

pub fn expand(input: TokenStream) -> TokenStream {
    match parse2::<MolInput>(input).and_then(codegen) {
        Ok(tokens) => tokens,
        Err(error) => error.to_compile_error(),
    }
}

fn codegen(input: MolInput) -> Result<TokenStream> {
    for path in &input.paths {
        for atom in path.atoms() {
            if let Atom::Port { name } = atom {
                return Err(Error::new(
                    name.span(),
                    "ports (`^name`) are only allowed in `frag!`",
                ));
            }
        }
    }

    let (specs, positions) = resolve_positions(&input.paths)?;

    let atoms = if specs.is_empty() {
        quote! {}
    } else {
        quote! { + atoms([ #(#specs),* ]) }
    };

    let mut bonds = Vec::new();
    for (path, row) in input.paths.iter().zip(&positions) {
        for (index, (op, _)) in path.rest.iter().enumerate() {
            let first = row[index].expect("mol! rejects ports before resolving");
            let second = row[index + 1].expect("mol! rejects ports before resolving");
            bonds.push(bond_term(op, first, second));
        }
    }

    Ok(quote! {
        {
            #[allow(unused_imports)]
            use ::umol_ast::ast::spec::{atoms, bond, double, single, triple, MoleculeSpec};
            (MoleculeSpec::new() #atoms #(+ #bonds)*).build()
        }
    })
}
