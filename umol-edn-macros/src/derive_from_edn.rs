//! `#[derive(FromEdn)]` proc macro for named structs and enums.
//!
//! ## Named structs
//!
//! Generates an `impl FromEdn<'de>` providing both `from_edn` (tree walk) and
//! `from_edn_str` (parser-deserializer fusion via `EdnStreamDeserializer`).
//! Default field rules:
//!
//! - `Option<T>`     → optional, missing key yields `None`.
//! - `Vec<T>`        → optional, missing key yields `Vec::new()`.
//! - `HashMap<K, V>` → optional, missing key yields `HashMap::new()`.
//! - `HashSet<T>`    → optional, missing key yields `HashSet::new()`.
//! - `#[edn(default)]` → optional, missing key yields `Default::default()`.
//! - everything else → required, missing key is an error.
//!
//! Field key naming: Rust `snake_case` is converted to EDN `kebab-case`
//! by default. Override with `#[edn(rename = "...")]`.
//!
//! ## Enums
//!
//! Generates an `impl FromEdn<'de>` that deserializes:
//! - Unit variants from EDN keywords: `:variant-name`
//! - Newtype variants from single-key maps: `{:variant-name value}`
//! - Tuple variants from single-key maps with vector values: `{:variant-name [v1 v2]}`
//! - Struct variants from single-key maps with map values: `{:variant-name {:field v}}`
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
            "FromEdn can only be derived on structs with named fields or enums",
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
    let mut bindings = Vec::new();
    let mut stream_decls = Vec::new();
    let mut stream_arms = Vec::new();
    let mut stream_finals = Vec::new();
    let mut field_idents = Vec::new();

    for field in fields {
        let field_info = parse_field(field)?;
        let ident = field_info.ident;
        let key = field_info.key;
        let kind = field_info.kind;

        let binding = match &kind {
            FieldKind::Option(inner) => quote! {
                let #ident: ::std::option::Option<#inner> = h.optional(#key)?;
            },
            FieldKind::Defaulted(ty) => quote! {
                let #ident: #ty = h
                    .optional::<#ty>(#key)?
                    .unwrap_or_default();
            },
            FieldKind::Required(ty) => quote! {
                let #ident: #ty = h.required(#key)?;
            },
        };

        let missing_key = key.clone();
        let (decl, arm, final_bind) = match &kind {
            FieldKind::Option(inner) => (
                quote! {
                    let mut #ident: ::std::option::Option<#inner> = ::std::option::Option::None;
                },
                quote! {
                    #key => {
                        let __slice = __de.read_value_slice()?;
                        #ident = ::std::option::Option::Some(
                            <#inner as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?,
                        );
                    }
                },
                quote! {},
            ),
            FieldKind::Defaulted(ty) => (
                quote! {
                    let mut #ident: #ty = <#ty as ::std::default::Default>::default();
                },
                quote! {
                    #key => {
                        let __slice = __de.read_value_slice()?;
                        #ident = <#ty as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?;
                    }
                },
                quote! {},
            ),
            FieldKind::Required(ty) => (
                quote! {
                    let mut #ident: ::std::option::Option<#ty> = ::std::option::Option::None;
                },
                quote! {
                    #key => {
                        let __slice = __de.read_value_slice()?;
                        #ident = ::std::option::Option::Some(
                            <#ty as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?,
                        );
                    }
                },
                quote! {
                    let #ident = #ident.ok_or_else(|| ::umol_edn::DeError::MissingField {
                        key: #missing_key.to_string(),
                        path: ::std::vec::Vec::new(),
                    })?;
                },
            ),
        };

        field_idents.push(ident);
        bindings.push(binding);
        stream_decls.push(decl);
        stream_arms.push(arm);
        stream_finals.push(final_bind);
    }

    let expected_label = format!("{} map", struct_name);

    Ok(quote! {
        impl<'de> ::umol_edn::FromEdn<'de> for #struct_name {
            fn from_edn(
                edn: &::umol_edn::Edn<'de>,
            ) -> ::std::result::Result<Self, ::umol_edn::DeError> {
                let m = match edn {
                    ::umol_edn::Edn::Map(m) => m,
                    other => {
                        return ::std::result::Result::Err(
                            ::umol_edn::DeError::TypeMismatch {
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

            fn from_edn_str(
                input: &'de str,
            ) -> ::std::result::Result<Self, ::umol_edn::EdnError> {
                let mut __de = ::umol_edn::EdnStreamDeserializer::new(input);
                #(#stream_decls)*
                __de.consume_byte(b'{')?;
                loop {
                    if __de.try_consume_byte(b'}')? {
                        break;
                    }
                    let __key = __de.read_keyword_name()?;
                    match __key.as_ref() {
                        #(#stream_arms)*
                        _ => __de.read_skip_value()?,
                    }
                }
                __de.expect_eof()?;
                #(#stream_finals)*
                ::std::result::Result::Ok(Self {
                    #(#field_idents),*
                })
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Enums
// ---------------------------------------------------------------------------

fn expand_enum(name: &syn::Ident, data: &DataEnum) -> Result<TokenStream2, syn::Error> {
    let is_all_unit = data
        .variants
        .iter()
        .all(|v| matches!(v.fields, Fields::Unit));

    let mut infos = Vec::new();
    for variant in &data.variants {
        let rename = read_rename(&variant.attrs)?;
        let key = rename.unwrap_or_else(|| variant.ident.to_string().to_kebab_case());
        infos.push(VariantInfo {
            ident: &variant.ident,
            key,
            fields: &variant.fields,
        });
    }

    if is_all_unit {
        expand_unit_enum(name, &infos)
    } else {
        expand_mixed_enum(name, &infos)
    }
}

struct VariantInfo<'a> {
    ident: &'a syn::Ident,
    key: String,
    fields: &'a Fields,
}

fn expand_unit_enum(
    name: &syn::Ident,
    variants: &[VariantInfo],
) -> Result<TokenStream2, syn::Error> {
    let from_arms: Vec<_> = variants
        .iter()
        .map(|v| {
            let ident = v.ident;
            let key = &v.key;
            quote! { #key => ::std::result::Result::Ok(#name::#ident), }
        })
        .collect();

    let stream_arms = from_arms.clone();
    let type_name = name.to_string();

    Ok(quote! {
        impl<'de> ::umol_edn::FromEdn<'de> for #name {
            fn from_edn(
                edn: &::umol_edn::Edn<'de>,
            ) -> ::std::result::Result<Self, ::umol_edn::DeError> {
                let s = match edn {
                    ::umol_edn::Edn::Keyword(k) => k.as_str(),
                    ::umol_edn::Edn::Str(s) => s.as_ref(),
                    other => return ::std::result::Result::Err(
                        ::umol_edn::DeError::TypeMismatch {
                            expected: "keyword",
                            got: other.kind(),
                            path: ::std::vec::Vec::new(),
                        },
                    ),
                };
                match s {
                    #(#from_arms)*
                    _ => ::std::result::Result::Err(::umol_edn::DeError::Custom(
                        format!("unknown {} variant: {:?}", #type_name, s),
                    )),
                }
            }

            fn from_edn_str(
                input: &'de str,
            ) -> ::std::result::Result<Self, ::umol_edn::EdnError> {
                let mut __de = ::umol_edn::EdnStreamDeserializer::new(input);
                let __kw = __de.read_keyword_name()?;
                __de.expect_eof()?;
                match __kw.as_ref() {
                    #(#stream_arms)*
                    _ => ::std::result::Result::Err(::umol_edn::EdnError::De(
                        ::umol_edn::DeError::Custom(
                            format!("unknown {} variant: {:?}", #type_name, __kw),
                        ),
                    )),
                }
            }
        }
    })
}

fn expand_mixed_enum(
    name: &syn::Ident,
    variants: &[VariantInfo],
) -> Result<TokenStream2, syn::Error> {
    let type_name = name.to_string();

    let mut unit_arms = Vec::new();
    let mut map_arms = Vec::new();

    for v in variants {
        let ident = v.ident;
        let key = &v.key;
        match v.fields {
            Fields::Unit => {
                unit_arms.push(quote! {
                    #key => ::std::result::Result::Ok(#name::#ident),
                });
                // Also accept unit variants in map form: {:variant nil}
                map_arms.push(quote! {
                    #key => ::std::result::Result::Ok(#name::#ident),
                });
            }
            Fields::Unnamed(fields) if fields.unnamed.len() == 1 => {
                let inner_ty = &fields.unnamed[0].ty;
                map_arms.push(quote! {
                    #key => {
                        let __inner = <#inner_ty as ::umol_edn::FromEdn<'de>>::from_edn(__val)?;
                        ::std::result::Result::Ok(#name::#ident(__inner))
                    }
                });
            }
            Fields::Unnamed(fields) => {
                let count = fields.unnamed.len();
                let field_parsers: Vec<_> = fields
                    .unnamed
                    .iter()
                    .enumerate()
                    .map(|(i, f)| {
                        let ty = &f.ty;
                        let idx = syn::Index::from(i);
                        quote! {
                            <#ty as ::umol_edn::FromEdn<'de>>::from_edn(&__seq[#idx])?
                        }
                    })
                    .collect();

                map_arms.push(quote! {
                    #key => {
                        let __seq = match __val {
                            ::umol_edn::Edn::Vector(s) => s,
                            other => return ::std::result::Result::Err(
                                ::umol_edn::DeError::TypeMismatch {
                                    expected: "vector",
                                    got: other.kind(),
                                    path: ::std::vec::Vec::new(),
                                },
                            ),
                        };
                        if __seq.len() != #count {
                            return ::std::result::Result::Err(::umol_edn::DeError::Custom(
                                format!(
                                    "expected {} elements for {} variant {}, got {}",
                                    #count, #type_name, #key, __seq.len(),
                                ),
                            ));
                        }
                        ::std::result::Result::Ok(#name::#ident(#(#field_parsers),*))
                    }
                });
            }
            Fields::Named(fields) => {
                let field_bindings: Vec<_> = fields
                    .named
                    .iter()
                    .map(|f| {
                        let fi = parse_field(f).unwrap();
                        let f_ident = fi.ident;
                        let f_key = fi.key;
                        match fi.kind {
                            FieldKind::Option(inner) => quote! {
                                let #f_ident: ::std::option::Option<#inner> =
                                    __h.optional(#f_key)?;
                            },
                            FieldKind::Defaulted(ty) => quote! {
                                let #f_ident: #ty =
                                    __h.optional::<#ty>(#f_key)?.unwrap_or_default();
                            },
                            FieldKind::Required(ty) => quote! {
                                let #f_ident: #ty = __h.required(#f_key)?;
                            },
                        }
                    })
                    .collect();
                let f_idents: Vec<_> = fields
                    .named
                    .iter()
                    .map(|f| f.ident.as_ref().unwrap())
                    .collect();

                map_arms.push(quote! {
                    #key => {
                        let __inner_map = match __val {
                            ::umol_edn::Edn::Map(m) => m,
                            other => return ::std::result::Result::Err(
                                ::umol_edn::DeError::TypeMismatch {
                                    expected: "map",
                                    got: other.kind(),
                                    path: ::std::vec::Vec::new(),
                                },
                            ),
                        };
                        let mut __h = ::umol_edn::EdnMapHelper::new(__inner_map);
                        #(#field_bindings)*
                        ::std::result::Result::Ok(#name::#ident { #(#f_idents),* })
                    }
                });
            }
        }
    }

    let has_unit_variants = !unit_arms.is_empty();
    let has_data_variants = map_arms.iter().any(|_| true);

    // Build the main match expression
    let keyword_branch = if has_unit_variants {
        quote! {
            ::umol_edn::Edn::Keyword(k) => {
                match k.as_str() {
                    #(#unit_arms)*
                    _ => ::std::result::Result::Err(::umol_edn::DeError::Custom(
                        format!("unknown {} variant: {:?}", #type_name, k.as_str()),
                    )),
                }
            }
        }
    } else {
        quote! {}
    };

    let map_branch = if has_data_variants {
        quote! {
            ::umol_edn::Edn::Map(m) => {
                if m.len() != 1 {
                    return ::std::result::Result::Err(::umol_edn::DeError::Custom(
                        format!("expected single-key map for {} variant, got {} keys",
                                #type_name, m.len()),
                    ));
                }
                let (__k, __val) = m.iter().next().unwrap();
                let __tag = match __k {
                    ::umol_edn::Edn::Keyword(k) => k.as_str(),
                    other => return ::std::result::Result::Err(
                        ::umol_edn::DeError::TypeMismatch {
                            expected: "keyword",
                            got: other.kind(),
                            path: ::std::vec::Vec::new(),
                        },
                    ),
                };
                match __tag {
                    #(#map_arms)*
                    _ => ::std::result::Result::Err(::umol_edn::DeError::Custom(
                        format!("unknown {} variant: {:?}", #type_name, __tag),
                    )),
                }
            }
        }
    } else {
        quote! {}
    };

    let expected = if has_unit_variants && has_data_variants {
        "keyword or single-key map"
    } else if has_unit_variants {
        "keyword"
    } else {
        "single-key map"
    };

    Ok(quote! {
        impl<'de> ::umol_edn::FromEdn<'de> for #name {
            fn from_edn(
                edn: &::umol_edn::Edn<'de>,
            ) -> ::std::result::Result<Self, ::umol_edn::DeError> {
                match edn {
                    #keyword_branch
                    #map_branch
                    other => ::std::result::Result::Err(
                        ::umol_edn::DeError::TypeMismatch {
                            expected: #expected,
                            got: other.kind(),
                            path: ::std::vec::Vec::new(),
                        },
                    ),
                }
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Field parsing helpers
// ---------------------------------------------------------------------------

struct FieldInfo {
    ident: syn::Ident,
    key: String,
    kind: FieldKind,
}

enum FieldKind {
    Option(Type),
    Defaulted(Type),
    Required(Type),
}

fn parse_field(field: &Field) -> Result<FieldInfo, syn::Error> {
    let ident = field
        .ident
        .clone()
        .ok_or_else(|| syn::Error::new_spanned(field, "field must be named"))?;

    let key = read_rename(&field.attrs)?.unwrap_or_else(|| to_kebab_case(&ident.to_string()));

    let has_default = has_default_attr(&field.attrs);

    const DEFAULTED_TYPES: &[&str] = &[
        "Vec", "HashMap", "HashSet", "BTreeMap", "BTreeSet", "IndexMap",
    ];

    let kind = if let Some(inner) = inner_of(&field.ty, "Option") {
        FieldKind::Option(inner.clone())
    } else if has_default || DEFAULTED_TYPES.iter().any(|t| is_type(&field.ty, t)) {
        FieldKind::Defaulted(field.ty.clone())
    } else {
        FieldKind::Required(field.ty.clone())
    };

    Ok(FieldInfo { ident, key, kind })
}

/// Read `#[edn(rename = "key")]` from attributes.
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
                Ok(()) // handled by has_default_attr
            } else {
                Err(meta.error("unsupported #[edn(...)] attribute"))
            }
        })?;
    }
    Ok(found)
}

/// Check for `#[edn(default)]` on a field.
fn has_default_attr(attrs: &[syn::Attribute]) -> bool {
    for attr in attrs {
        if !attr.path().is_ident("edn") {
            continue;
        }
        let mut found = false;
        let _ = attr.parse_nested_meta(|meta| {
            if meta.path.is_ident("default") {
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

/// Check if `ty` is `Name<...>`.
fn is_type(ty: &Type, name: &str) -> bool {
    if let Type::Path(type_path) = ty {
        if let Some(segment) = type_path.path.segments.last() {
            return segment.ident == name;
        }
    }
    false
}

/// Convert `snake_case` to `kebab-case`.
fn to_kebab_case(s: &str) -> String {
    s.replace('_', "-")
}

