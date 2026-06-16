//! Proc macros for `umol-ast`. Exports `#[derive(Lattice)]`, which generates a
//! field-wise `Lattice` impl for structs whose fields all themselves implement
//! `Lattice`.

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, Data, DataStruct, DeriveInput, Fields, FieldsNamed};

/// Derive `Lattice` for a struct by propagating the trait method-by-method
/// across the named fields. Every field type must itself implement `Lattice`
/// (the impl emits no extra bounds; missing impls surface as compile errors at
/// the call site).
///
/// Generated impl:
/// - `is_undetermined` / `is_ground` / `matches`: conjunction of the field-wise calls
/// - `meet`: field-wise meet; returns `None` as soon as any field's meet is `None`
/// - `join`: field-wise join
///
/// Unsupported: tuple structs, unit structs, enums, generic types. These produce
/// a compile-time error.
#[proc_macro_derive(Lattice)]
pub fn derive_lattice(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;
    let generics = &input.generics;
    if !generics.params.is_empty() {
        return syn::Error::new_spanned(
            generics,
            "derive(Lattice) does not support generic types yet",
        )
        .into_compile_error()
        .into();
    }
    let fields = match &input.data {
        Data::Struct(DataStruct {
            fields: Fields::Named(FieldsNamed { named, .. }),
            ..
        }) => named,
        Data::Struct(DataStruct { fields, .. }) => {
            return syn::Error::new_spanned(
                fields,
                "derive(Lattice) requires a struct with named fields",
            )
            .into_compile_error()
            .into();
        }
        Data::Enum(_) | Data::Union(_) => {
            return syn::Error::new_spanned(
                name,
                "derive(Lattice) only supports structs; enum/union impls must be hand-rolled",
            )
            .into_compile_error()
            .into();
        }
    };
    let field_names: Vec<&syn::Ident> = fields
        .iter()
        .map(|f| f.ident.as_ref().expect("named field"))
        .collect();
    let lattice = quote!(crate::ast::Lattice);
    quote! {
        impl #lattice for #name {
            fn is_undetermined(&self) -> bool {
                true #( && #lattice::is_undetermined(&self.#field_names) )*
            }
            fn is_ground(&self) -> bool {
                true #( && #lattice::is_ground(&self.#field_names) )*
            }
            fn meet(&self, other: &Self) -> ::core::option::Option<Self> {
                ::core::option::Option::Some(Self {
                    #( #field_names: #lattice::meet(&self.#field_names, &other.#field_names)?, )*
                })
            }
            fn join(&self, other: &Self) -> Self {
                Self {
                    #( #field_names: #lattice::join(&self.#field_names, &other.#field_names), )*
                }
            }
            fn matches(&self, target: &Self) -> bool {
                true #( && #lattice::matches(&self.#field_names, &target.#field_names) )*
            }
        }
    }
    .into()
}
