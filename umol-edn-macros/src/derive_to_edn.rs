//! `#[derive(ToEdn)]` proc macro for named structs and enums.
//!
//! ## Named structs
//!
//! Generates an `impl ToEdn` that builds an `EdnMap` from each field, using
//! kebab-case field keys by default. `Option<T>` fields are skipped when
//! `None`; all other fields (including empty `Vec<T>`) are emitted.
//!
//! ## Enums
//!
//! Generates an `impl ToEdn` that serializes:
//! - Unit variants as EDN keywords: `:variant-name`
//! - Newtype variants as single-key maps: `{:variant-name value}`
//! - Tuple variants as single-key maps with vector values: `{:variant-name [v1 v2]}`
//! - Struct variants as single-key maps with map values: `{:variant-name {:field v}}`
//!
//! Variant naming: Rust `PascalCase` is converted to `kebab-case` by default.
//! Override with `#[edn(rename = "...")]`.

use heck::ToKebabCase;
use proc_macro2::TokenStream as TokenStream2;
use quote::quote;
use syn::{Data, DataEnum, DataStruct, DeriveInput, Field, Fields, GenericArgument, PathArguments, Type};

pub fn expand(input: DeriveInput) -> Result<TokenStream2, syn::Error> {
    match &input.data {
        Data::Struct(DataStruct {
            fields: Fields::Named(named),
            ..
        }) => expand_struct(&input.ident, &named.named),
        Data::Enum(data_enum) => expand_enum(&input.ident, data_enum),
        _ => Err(syn::Error::new_spanned(
            input.ident,
            "ToEdn can only be derived on structs with named fields or enums",
        )),
    }
}

// ---------------------------------------------------------------------------
// Structs
// ---------------------------------------------------------------------------

fn expand_struct(
    struct_name: &syn::Ident,
    fields: &syn::punctuated::Punctuated<Field, syn::Token![,]>,
) -> Result<TokenStream2, syn::Error> {
    let mut inserts = Vec::new();
    let len = fields.len();

    for field in fields {
        let info = parse_field(field)?;
        let ident = info.ident;
        let key = info.key;

        let insert = if info.is_option {
            quote! {
                if let ::std::option::Option::Some(__v) = &self.#ident {
                    m.insert(
                        ::umol_edn::Edn::keyword(#key),
                        ::umol_edn::ToEdn::to_edn(__v),
                    );
                }
            }
        } else {
            quote! {
                m.insert(
                    ::umol_edn::Edn::keyword(#key),
                    ::umol_edn::ToEdn::to_edn(&self.#ident),
                );
            }
        };

        inserts.push(insert);
    }

    Ok(quote! {
        impl ::umol_edn::ToEdn for #struct_name {
            fn to_edn(&self) -> ::umol_edn::Edn<'static> {
                let mut m = ::umol_edn::EdnMap::with_capacity(#len);
                #(#inserts)*
                ::umol_edn::Edn::Map(m)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

fn expand_enum(name: &syn::Ident, data: &DataEnum) -> Result<TokenStream2, syn::Error> {
    let mut arms = Vec::new();

    for variant in &data.variants {
        let ident = &variant.ident;
        let rename = read_rename(&variant.attrs)?;
        let key = rename.unwrap_or_else(|| ident.to_string().to_kebab_case());

        let arm = match &variant.fields {
            Fields::Unit => {
                quote! {
                    #name::#ident => ::umol_edn::Edn::keyword(#key),
                }
            }
            Fields::Unnamed(fields) if fields.unnamed.len() == 1 => {
                quote! {
                    #name::#ident(ref __v0) => {
                        let mut __m = ::umol_edn::EdnMap::with_capacity(1);
                        __m.insert(
                            ::umol_edn::Edn::keyword(#key),
                            ::umol_edn::ToEdn::to_edn(__v0),
                        );
                        ::umol_edn::Edn::Map(__m)
                    }
                }
            }
            Fields::Unnamed(fields) => {
                let field_refs: Vec<_> = (0..fields.unnamed.len())
                    .map(|i| {
                        let name = syn::Ident::new(&format!("__v{}", i), ident.span());
                        quote! { ref #name }
                    })
                    .collect();
                let field_edns: Vec<_> = (0..fields.unnamed.len())
                    .map(|i| {
                        let name = syn::Ident::new(&format!("__v{}", i), ident.span());
                        quote! { ::umol_edn::ToEdn::to_edn(#name) }
                    })
                    .collect();

                quote! {
                    #name::#ident(#(#field_refs),*) => {
                        let mut __m = ::umol_edn::EdnMap::with_capacity(1);
                        __m.insert(
                            ::umol_edn::Edn::keyword(#key),
                            ::umol_edn::Edn::Vector(
                                vec![#(#field_edns),*].into(),
                            ),
                        );
                        ::umol_edn::Edn::Map(__m)
                    }
                }
            }
            Fields::Named(fields) => {
                let field_refs: Vec<_> = fields
                    .named
                    .iter()
                    .map(|f| {
                        let f_ident = f.ident.as_ref().unwrap();
                        quote! { ref #f_ident }
                    })
                    .collect();

                let field_count = fields.named.len();
                let field_inserts: Vec<_> = fields
                    .named
                    .iter()
                    .map(|f| {
                        let f_ident = f.ident.as_ref().unwrap();
                        let f_key = read_rename(&f.attrs)
                            .unwrap()
                            .unwrap_or_else(|| to_kebab_case(&f_ident.to_string()));
                        let is_opt = is_option_type(&f.ty);

                        if is_opt {
                            quote! {
                                if let ::std::option::Option::Some(__v) = #f_ident {
                                    __inner.insert(
                                        ::umol_edn::Edn::keyword(#f_key),
                                        ::umol_edn::ToEdn::to_edn(__v),
                                    );
                                }
                            }
                        } else {
                            quote! {
                                __inner.insert(
                                    ::umol_edn::Edn::keyword(#f_key),
                                    ::umol_edn::ToEdn::to_edn(#f_ident),
                                );
                            }
                        }
                    })
                    .collect();

                quote! {
                    #name::#ident { #(#field_refs),* } => {
                        let mut __inner = ::umol_edn::EdnMap::with_capacity(#field_count);
                        #(#field_inserts)*
                        let mut __m = ::umol_edn::EdnMap::with_capacity(1);
                        __m.insert(
                            ::umol_edn::Edn::keyword(#key),
                            ::umol_edn::Edn::Map(__inner),
                        );
                        ::umol_edn::Edn::Map(__m)
                    }
                }
            }
        };

        arms.push(arm);
    }

    Ok(quote! {
        impl ::umol_edn::ToEdn for #name {
            fn to_edn(&self) -> ::umol_edn::Edn<'static> {
                match self {
                    #(#arms)*
                }
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

struct FieldInfo {
    ident: syn::Ident,
    key: String,
    is_option: bool,
}

fn parse_field(field: &Field) -> Result<FieldInfo, syn::Error> {
    let ident = field
        .ident
        .clone()
        .ok_or_else(|| syn::Error::new_spanned(field, "field must be named"))?;

    let key = read_rename(&field.attrs)?.unwrap_or_else(|| to_kebab_case(&ident.to_string()));
    let is_option = is_option_type(&field.ty);

    Ok(FieldInfo {
        ident,
        key,
        is_option,
    })
}

fn read_rename(attrs: &[syn::Attribute]) -> Result<Option<String>, syn::Error> {
    let mut found: Option<String> = None;
    for attr in attrs {
        if !attr.path().is_ident("edn") {
            continue;
        }
        attr.parse_nested_meta(|meta| {
            if meta.path.is_ident("rename") {
                let value: syn::LitStr = meta.value()?.parse()?;
                found = Some(value.value());
                Ok(())
            } else if meta.path.is_ident("default") {
                Ok(()) // only relevant for FromEdn
            } else {
                Err(meta.error("unsupported #[edn(...)] attribute"))
            }
        })?;
    }
    Ok(found)
}

fn is_option_type(ty: &Type) -> bool {
    if let Type::Path(type_path) = ty {
        if let Some(segment) = type_path.path.segments.last() {
            if segment.ident == "Option" {
                if let PathArguments::AngleBracketed(args) = &segment.arguments {
                    if matches!(args.args.first(), Some(GenericArgument::Type(_))) {
                        return true;
                    }
                }
            }
        }
    }
    false
}

fn to_kebab_case(s: &str) -> String {
    s.replace('_', "-")
}

