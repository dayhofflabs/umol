//! `#[derive(ToEdn)]` proc macro for structs and enums.
//!
//! ## Structs
//!
//! Generates an `impl ToEdn` that builds an `EdnMap` from each field.
//!
//! Field key naming: Rust `snake_case` → EDN `kebab-case` by default.
//!
//! ## Container attributes
//!
//! - `#[edn(transparent)]` — single-field struct (named or tuple) delegates
//!   to the inner type's `ToEdn`.
//!
//! ## Field attributes
//!
//! - `#[edn(rename = "key")]` — override the EDN key for this field.
//! - `#[edn(skip)]` — exclude from serialization entirely.
//! - `#[edn(skip_if = "path::to::fn")]` — omit when the predicate returns
//!   true. The predicate receives `&FieldType`.
//!
//! Serialization rules (first match wins):
//!
//! - `#[edn(skip)]`      → never emitted.
//! - `#[edn(skip_if)]`   → emitted unless the predicate returns true.
//! - `Option<T>`          → emitted when `Some`, omitted when `None`.
//! - everything else      → always emitted.
//!
//! ## Enums
//!
//! Generates an `impl ToEdn` that serializes:
//! - Unit variants as EDN keywords: `:variant-name`
//! - Newtype variants as single-key maps: `{:variant-name value}`
//! - Tuple variants as single-key maps with vector values: `{:variant-name [v1 v2]}`
//! - Struct variants as single-key maps with map values: `{:variant-name {:field v}}`
//!
//! Variant naming: Rust `PascalCase` → `kebab-case` by default.
//! Override with `#[edn(rename = "...")]`.

use heck::ToKebabCase;
use proc_macro2::TokenStream as TokenStream2;
use quote::quote;
use syn::{Data, DataEnum, DataStruct, DeriveInput, Field, Fields, GenericArgument, PathArguments, Type};

pub fn expand(input: DeriveInput) -> Result<TokenStream2, syn::Error> {
    if has_container_attr(&input.attrs, "transparent") {
        return expand_transparent(&input);
    }

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

// Structs
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

        let insert = match &info.ser {
            FieldSer::Skip => continue,
            FieldSer::Option => quote! {
                if let ::std::option::Option::Some(__v) = &self.#ident {
                    m.insert(
                        ::umol_edn::Edn::keyword(#key),
                        ::umol_edn::ToEdn::to_edn(__v),
                    );
                }
            },
            FieldSer::SkipIf(predicate) => quote! {
                if !#predicate(&self.#ident) {
                    m.insert(
                        ::umol_edn::Edn::keyword(#key),
                        ::umol_edn::ToEdn::to_edn(&self.#ident),
                    );
                }
            },
            FieldSer::Normal => quote! {
                m.insert(
                    ::umol_edn::Edn::keyword(#key),
                    ::umol_edn::ToEdn::to_edn(&self.#ident),
                );
            },
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

// Enums
fn expand_enum(name: &syn::Ident, data: &DataEnum) -> Result<TokenStream2, syn::Error> {
    let mut arms = Vec::new();

    for variant in &data.variants {
        let ident = &variant.ident;
        let rename = FieldAttrs::parse(&variant.attrs)?.rename;
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
                let parsed_fields: Vec<_> = fields
                    .named
                    .iter()
                    .map(parse_field)
                    .collect::<Result<Vec<_>, _>>()?;

                let mut field_refs = Vec::new();
                let mut field_inserts = Vec::new();
                let field_count = parsed_fields.len();

                for fi in &parsed_fields {
                    let f_ident = &fi.ident;
                    let f_key = &fi.key;

                    match &fi.ser {
                        FieldSer::Skip => {
                            field_refs.push(quote! { #f_ident: _ });
                            continue;
                        }
                        FieldSer::Option => {
                            field_refs.push(quote! { ref #f_ident });
                            field_inserts.push(quote! {
                                if let ::std::option::Option::Some(__v) = #f_ident {
                                    __inner.insert(
                                        ::umol_edn::Edn::keyword(#f_key),
                                        ::umol_edn::ToEdn::to_edn(__v),
                                    );
                                }
                            });
                        }
                        FieldSer::SkipIf(predicate) => {
                            field_refs.push(quote! { ref #f_ident });
                            field_inserts.push(quote! {
                                if !#predicate(#f_ident) {
                                    __inner.insert(
                                        ::umol_edn::Edn::keyword(#f_key),
                                        ::umol_edn::ToEdn::to_edn(#f_ident),
                                    );
                                }
                            });
                        }
                        FieldSer::Normal => {
                            field_refs.push(quote! { ref #f_ident });
                            field_inserts.push(quote! {
                                __inner.insert(
                                    ::umol_edn::Edn::keyword(#f_key),
                                    ::umol_edn::ToEdn::to_edn(#f_ident),
                                );
                            });
                        }
                    }
                }

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

// Transparent structs
fn expand_transparent(input: &DeriveInput) -> Result<TokenStream2, syn::Error> {
    let name = &input.ident;
    let accessor = match &input.data {
        Data::Struct(DataStruct { fields: Fields::Named(named), .. }) => {
            if named.named.len() != 1 {
                return Err(syn::Error::new_spanned(
                    name,
                    "transparent requires exactly one field",
                ));
            }
            let ident = named.named.first().unwrap().ident.as_ref().unwrap();
            quote! { &self.#ident }
        }
        Data::Struct(DataStruct { fields: Fields::Unnamed(unnamed), .. }) => {
            if unnamed.unnamed.len() != 1 {
                return Err(syn::Error::new_spanned(
                    name,
                    "transparent requires exactly one field",
                ));
            }
            quote! { &self.0 }
        }
        _ => {
            return Err(syn::Error::new_spanned(
                name,
                "transparent can only be used on single-field structs",
            ));
        }
    };

    Ok(quote! {
        impl ::umol_edn::ToEdn for #name {
            fn to_edn(&self) -> ::umol_edn::Edn<'static> {
                ::umol_edn::ToEdn::to_edn(#accessor)
            }
        }
    })
}

fn has_container_attr(attrs: &[syn::Attribute], name: &str) -> bool {
    for attr in attrs {
        if !attr.path().is_ident("edn") {
            continue;
        }
        let mut found = false;
        let _ = attr.parse_nested_meta(|meta| {
            if meta.path.is_ident(name) {
                found = true;
            }
            Ok(())
        });
        if found {
            return true;
        }
    }
    false
}

// Helpers
enum FieldSer {
    Normal,
    Option,
    Skip,
    SkipIf(syn::Path),
}

struct FieldInfo {
    ident: syn::Ident,
    key: String,
    ser: FieldSer,
}

fn parse_field(field: &Field) -> Result<FieldInfo, syn::Error> {
    let ident = field
        .ident
        .clone()
        .ok_or_else(|| syn::Error::new_spanned(field, "field must be named"))?;

    let attrs = FieldAttrs::parse(&field.attrs)?;
    let key = attrs.rename.unwrap_or_else(|| ident.to_string().to_kebab_case());

    let ser = if attrs.skip {
        FieldSer::Skip
    } else if let Some(path) = attrs.skip_if {
        FieldSer::SkipIf(path)
    } else if is_option_type(&field.ty) {
        FieldSer::Option
    } else {
        FieldSer::Normal
    };

    Ok(FieldInfo { ident, key, ser })
}

struct FieldAttrs {
    rename: Option<String>,
    skip: bool,
    skip_if: Option<syn::Path>,
}

impl FieldAttrs {
    fn parse(attrs: &[syn::Attribute]) -> Result<Self, syn::Error> {
        let mut rename = None;
        let mut skip = false;
        let mut skip_if = None;
        for attr in attrs {
            if !attr.path().is_ident("edn") {
                continue;
            }
            attr.parse_nested_meta(|meta| {
                if meta.path.is_ident("rename") {
                    let value: syn::LitStr = meta.value()?.parse()?;
                    rename = Some(value.value());
                    Ok(())
                } else if meta.path.is_ident("default") {
                    Ok(()) // only relevant for FromEdn
                } else if meta.path.is_ident("skip") {
                    skip = true;
                    Ok(())
                } else if meta.path.is_ident("skip_if") {
                    let value: syn::LitStr = meta.value()?.parse()?;
                    skip_if = Some(value.parse()?);
                    Ok(())
                } else {
                    Err(meta.error("unsupported #[edn(...)] attribute"))
                }
            })?;
        }
        Ok(FieldAttrs { rename, skip, skip_if })
    }
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


