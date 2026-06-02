//! Proc macros for `umol-ast`. Currently exports `#[derive(Lattice)]`, which
//! generates a field-wise `Lattice` impl for structs whose fields all
//! themselves implement `Lattice`. An optional
//! `#[lattice(saturate = "fn_name")]` attribute wires a cross-field
//! propagation hook that runs at the end of every `meet`.

use proc_macro::TokenStream;
use quote::quote;
use syn::{
    parse_macro_input, Attribute, Data, DataStruct, DeriveInput, Expr, ExprLit, Fields,
    FieldsNamed, Lit,
};

/// Derive `Lattice` for a struct by propagating the trait method-by-method
/// across the named fields. Every field type must itself implement
/// `Lattice` (the impl emits no extra bounds; missing impls surface as
/// compile errors at the call site).
///
/// Generated impl:
/// - `is_undetermined` / `is_ground` / `matches`: conjunction of the field-wise calls
/// - `meet`: field-wise meet; returns `None` as soon as any field's meet is `None`; then calls `saturate` on the result and converts `Err(Contradiction)` to `None`
/// - `join`: field-wise join
/// - `saturate` is overridden iff the `#[lattice(saturate = "fn_name")]` attribute is present; otherwise the trait's no-op default applies
///
/// Unsupported: tuple structs, unit structs, enums, generic types. These
/// produce a compile-time error.
#[proc_macro_derive(Lattice, attributes(lattice))]
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
    let saturate_fn = match parse_saturate_attr(&input.attrs) {
        Ok(opt) => opt,
        Err(e) => return e.into_compile_error().into(),
    };
    let lattice = quote!(crate::ast::Lattice);
    let contradiction = quote!(crate::ast::Contradiction);
    let saturate_override = saturate_fn.map(|fn_ident| {
        quote! {
            fn saturate(&mut self) -> ::core::result::Result<(), #contradiction> {
                #fn_ident(self)
            }
        }
    });
    quote! {
        impl #lattice for #name {
            fn is_undetermined(&self) -> bool {
                true #( && #lattice::is_undetermined(&self.#field_names) )*
            }
            fn is_ground(&self) -> bool {
                true #( && #lattice::is_ground(&self.#field_names) )*
            }
            fn meet(&self, other: &Self) -> ::core::option::Option<Self> {
                let mut result = Self {
                    #( #field_names: #lattice::meet(&self.#field_names, &other.#field_names)?, )*
                };
                #lattice::saturate(&mut result).ok()?;
                ::core::option::Option::Some(result)
            }
            fn join(&self, other: &Self) -> Self {
                Self {
                    #( #field_names: #lattice::join(&self.#field_names, &other.#field_names), )*
                }
            }
            fn matches(&self, target: &Self) -> bool {
                true #( && #lattice::matches(&self.#field_names, &target.#field_names) )*
            }
            #saturate_override
        }
    }
    .into()
}

/// Parse `#[lattice(saturate = "fn_name")]` from the struct attributes.
/// Returns `Ok(Some(ident))` if present and well-formed, `Ok(None)` if
/// absent, `Err` on malformed attribute.
fn parse_saturate_attr(attrs: &[Attribute]) -> syn::Result<Option<syn::Ident>> {
    for attr in attrs {
        if !attr.path().is_ident("lattice") {
            continue;
        }
        let mut result: Option<syn::Ident> = None;
        attr.parse_nested_meta(|meta| {
            if meta.path.is_ident("saturate") {
                let value: Expr = meta.value()?.parse()?;
                if let Expr::Lit(ExprLit {
                    lit: Lit::Str(s), ..
                }) = value
                {
                    let ident: syn::Ident = syn::parse_str(&s.value())?;
                    result = Some(ident);
                    return Ok(());
                }
                Err(meta.error("expected `saturate = \"fn_name\"`"))
            } else {
                Err(meta.error("unknown lattice attribute key"))
            }
        })?;
        return Ok(result);
    }
    Ok(None)
}
