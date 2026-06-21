//! Proc macros for `umol-ast`: `#[derive(Lattice)]` and `#[derive(Canonicalize)]`,
//! each generated field-wise over a struct's named fields (every field type must
//! itself implement the trait).

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, Data, DataStruct, DeriveInput, Fields, FieldsNamed, Ident};

/// Named fields of a struct, or a compile-error `TokenStream` for tuple/unit
/// structs, enums, unions, and generic types. `derive` names the trait for the
/// error message.
fn named_struct_fields(input: &DeriveInput, derive: &str) -> Result<Vec<Ident>, TokenStream> {
    if !input.generics.params.is_empty() {
        return Err(syn::Error::new_spanned(
            &input.generics,
            format!("derive({derive}) does not support generic types yet"),
        )
        .into_compile_error()
        .into());
    }
    match &input.data {
        Data::Struct(DataStruct {
            fields: Fields::Named(FieldsNamed { named, .. }),
            ..
        }) => Ok(named
            .iter()
            .map(|f| f.ident.clone().expect("named field"))
            .collect()),
        Data::Struct(DataStruct { fields, .. }) => Err(syn::Error::new_spanned(
            fields,
            format!("derive({derive}) requires a struct with named fields"),
        )
        .into_compile_error()
        .into()),
        Data::Enum(_) | Data::Union(_) => Err(syn::Error::new_spanned(
            &input.ident,
            format!("derive({derive}) only supports structs; enum/union impls must be hand-rolled"),
        )
        .into_compile_error()
        .into()),
    }
}

/// Derive `Lattice` field-wise across the named fields.
///
/// - `is_undetermined` / `is_ground`: conjunction of the field-wise calls
/// - `meet`: field-wise meet; returns `None` as soon as any field's meet is `None`
/// - `join`: field-wise join
/// - `matches`: conjunction of the field-wise `matches`. Equal to the trait's
///   `meet`-derived default (a struct of canonical fields is canonical), but built
///   from each field's `matches` directly so the per-candidate path allocates no
///   intermediate `meet`/`canonical`.
#[proc_macro_derive(Lattice)]
pub fn derive_lattice(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let fields = match named_struct_fields(&input, "Lattice") {
        Ok(f) => f,
        Err(e) => return e,
    };
    let name = &input.ident;
    let lattice = quote!(crate::ast::Lattice);
    quote! {
        impl #lattice for #name {
            fn is_undetermined(&self) -> bool {
                true #( && #lattice::is_undetermined(&self.#fields) )*
            }
            fn is_ground(&self) -> bool {
                true #( && #lattice::is_ground(&self.#fields) )*
            }
            fn meet(&self, other: &Self) -> ::core::option::Option<Self> {
                ::core::option::Option::Some(Self {
                    #( #fields: #lattice::meet(&self.#fields, &other.#fields)?, )*
                })
            }
            fn join(&self, other: &Self) -> Self {
                Self {
                    #( #fields: #lattice::join(&self.#fields, &other.#fields), )*
                }
            }
            fn matches(&self, target: &Self) -> bool {
                true #( && #lattice::matches(&self.#fields, &target.#fields) )*
            }
        }
    }
    .into()
}

/// Derive `Canonicalize` by canonicalizing each named field. Every field type
/// must itself implement `Canonicalize`. `canonical` uses the trait default
/// (the fast-path borrow is a per-type override, not generated here).
#[proc_macro_derive(Canonicalize)]
pub fn derive_canonicalize(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let fields = match named_struct_fields(&input, "Canonicalize") {
        Ok(f) => f,
        Err(e) => return e,
    };
    let name = &input.ident;
    let canon = quote!(crate::ast::Canonicalize);
    let contradiction = quote!(crate::ast::Contradiction);
    quote! {
        impl #canon for #name {
            fn canonicalize(self) -> ::core::result::Result<Self, #contradiction> {
                ::core::result::Result::Ok(Self {
                    #( #fields: #canon::canonicalize(self.#fields)?, )*
                })
            }
        }
    }
    .into()
}
