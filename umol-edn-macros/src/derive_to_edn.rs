//! `#[derive(ToEdn)]` proc macro for named structs.
//!
//! Generates an `impl ToEdn` that builds an `EdnMap` from each field, using
//! kebab-case field keys by default. `Option<T>` fields are skipped when
//! `None`; all other fields (including empty `Vec<T>`) are emitted.

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
                "ToEdn can only be derived on structs with named fields",
            ));
        }
    };

    let mut inserts = Vec::new();
    let len = fields.len();

    for field in &fields {
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
            fn to_edn(&self) -> ::umol_edn::Edn<'_> {
                let mut m = ::umol_edn::EdnMap::with_capacity(#len);
                #(#inserts)*
                ::umol_edn::Edn::Map(m)
            }
        }
    })
}

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

    let key = read_rename_attr(field)?.unwrap_or_else(|| to_kebab_case(&ident.to_string()));
    let is_option = is_option_type(&field.ty);

    Ok(FieldInfo {
        ident,
        key,
        is_option,
    })
}

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
