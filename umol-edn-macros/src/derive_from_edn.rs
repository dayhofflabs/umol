//! `#[derive(FromEdn)]` proc macro for named structs.
//!
//! Generates an `impl FromEdn<'de>` that walks an `Edn::Map`, calling
//! `EdnMapHelper::required` / `optional` for each field. The default
//! field rules are:
//!
//! - `Option<T>` → optional, missing key yields `None`.
//! - `Vec<T>`    → optional, missing key yields `Vec::new()`.
//! - everything else → required, missing key is an error.
//!
//! Field key naming: Rust `snake_case` is converted to EDN `kebab-case`
//! by default. Override with `#[edn(rename = "...")]`.

use proc_macro2::TokenStream as TokenStream2;
use quote::quote;
use syn::{Data, DataStruct, DeriveInput, Field, Fields, GenericArgument, PathArguments, Type};

pub fn expand(input: DeriveInput) -> Result<TokenStream2, syn::Error> {
    let struct_name = input.ident.clone();

    let fields = match input.data {
        Data::Struct(DataStruct {
            fields: Fields::Named(named),
            ..
        }) => named.named,
        _ => {
            return Err(syn::Error::new_spanned(
                struct_name,
                "FromEdn can only be derived on structs with named fields",
            ));
        }
    };

    let mut bindings = Vec::new();
    let mut field_idents = Vec::new();

    for field in &fields {
        let field_info = parse_field(field)?;
        let ident = field_info.ident;
        let key = field_info.key;
        let kind = field_info.kind;

        let binding = match kind {
            FieldKind::Option(inner) => quote! {
                let #ident: ::std::option::Option<#inner> = h.optional(#key)?;
            },
            FieldKind::Vec(inner) => quote! {
                let #ident: ::std::vec::Vec<#inner> = h
                    .optional::<::std::vec::Vec<#inner>>(#key)?
                    .unwrap_or_default();
            },
            FieldKind::Required(ty) => quote! {
                let #ident: #ty = h.required(#key)?;
            },
        };

        field_idents.push(ident);
        bindings.push(binding);
    }

    let expected_label = format!("{} map", struct_name);

    Ok(quote! {
        impl<'de> ::umol_edn::FromEdn<'de> for #struct_name {
            fn from_edn(
                edn: &::umol_edn::Edn<'de>,
            ) -> ::std::result::Result<Self, ::umol_edn::EdnError> {
                let m = match edn {
                    ::umol_edn::Edn::Map(m) => m,
                    other => {
                        return ::std::result::Result::Err(
                            ::umol_edn::EdnError::TypeMismatch {
                                expected: #expected_label,
                                got: other.kind(),
                                path: ::std::vec::Vec::new(),
                            },
                        );
                    }
                };
                let mut h = ::umol_edn::EdnMapHelper::new(m);
                #(#bindings)*
                ::std::result::Result::Ok(Self {
                    #(#field_idents),*
                })
            }
        }
    })
}

struct FieldInfo {
    ident: syn::Ident,
    key: String,
    kind: FieldKind,
}

enum FieldKind {
    Option(Type),
    Vec(Type),
    Required(Type),
}

fn parse_field(field: &Field) -> Result<FieldInfo, syn::Error> {
    let ident = field
        .ident
        .clone()
        .ok_or_else(|| syn::Error::new_spanned(field, "field must be named"))?;

    let key = read_rename_attr(field)?.unwrap_or_else(|| to_kebab_case(&ident.to_string()));

    let kind = if let Some(inner) = inner_of(&field.ty, "Option") {
        FieldKind::Option(inner.clone())
    } else if let Some(inner) = inner_of(&field.ty, "Vec") {
        FieldKind::Vec(inner.clone())
    } else {
        FieldKind::Required(field.ty.clone())
    };

    Ok(FieldInfo { ident, key, kind })
}

/// Read `#[edn(rename = "key")]` from a field's attributes, if present.
fn read_rename_attr(field: &Field) -> Result<Option<String>, syn::Error> {
    let mut found: Option<String> = None;
    for attr in &field.attrs {
        if !attr.path().is_ident("edn") {
            continue;
        }
        attr.parse_nested_meta(|meta| {
            if meta.path.is_ident("rename") {
                let value: syn::LitStr = meta.value()?.parse()?;
                found = Some(value.value());
                Ok(())
            } else {
                Err(meta.error("unsupported #[edn(...)] attribute"))
            }
        })?;
    }
    Ok(found)
}

/// If `ty` is `Wrapper<T>`, return `T`. Otherwise return `None`.
fn inner_of<'a>(ty: &'a Type, wrapper: &str) -> Option<&'a Type> {
    if let Type::Path(type_path) = ty {
        if let Some(segment) = type_path.path.segments.last() {
            if segment.ident == wrapper {
                if let PathArguments::AngleBracketed(args) = &segment.arguments {
                    if let Some(GenericArgument::Type(inner)) = args.args.first() {
                        return Some(inner);
                    }
                }
            }
        }
    }
    None
}

/// Convert `snake_case` to `kebab-case`.
fn to_kebab_case(s: &str) -> String {
    s.replace('_', "-")
}
