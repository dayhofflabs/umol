//! `#[derive(FromEdn)]` proc macro for structs and enums.
//!
//! ## Structs
//!
//! Generates an `impl FromEdn<'de>` providing both `from_edn` (tree walk) and
//! `from_edn_str` (parser-deserializer fusion via `EdnStreamDeserializer`).
//!
//! Field key naming: Rust `snake_case` → EDN `kebab-case` by default.
//!
//! ## Container attributes
//!
//! - `#[edn(transparent)]` — single-field struct (named or tuple) delegates
//!   to the inner type's `FromEdn`.
//! - `#[edn(deny_unknown_fields)]` — error on unrecognized map keys in both
//!   `from_edn` and `from_edn_str`.
//! - `#[edn(default)]` — all fields become optional, using
//!   `Default::default()` when missing.
//!
//! ## Field attributes
//!
//! - `#[edn(rename = "key")]` — override the EDN key for this field.
//! - `#[edn(default)]` — use `Default::default()` when this key is missing.
//! - `#[edn(skip)]` — exclude from deserialization; always `Default::default()`.
//!
//! Field rules (first match wins):
//!
//! - `#[edn(skip)]`     → always `Default::default()`, key ignored if present.
//! - `Option<T>`         → optional, missing key yields `None`.
//! - `#[edn(default)]` (field or container) → missing key yields `Default::default()`.
//! - everything else     → required, missing key is an error.
//!
//! ## Enums
//!
//! Generates an `impl FromEdn<'de>` with both `from_edn` and `from_edn_str`:
//! - Unit variants from EDN keywords: `:variant-name`
//! - Newtype variants from single-key maps: `{:variant-name value}`
//! - Tuple variants from single-key maps with vector values: `{:variant-name [v1 v2]}`
//! - Struct variants from single-key maps with map values: `{:variant-name {:field v}}`
//!
//! Variant naming: Rust `PascalCase` → `kebab-case` by default.
//! Override with `#[edn(rename = "...")]`.

use heck::ToKebabCase;
use proc_macro2::TokenStream as TokenStream2;
use quote::quote;
use syn::{Data, DataEnum, DataStruct, DeriveInput, Field, Fields, GenericArgument, PathArguments, Type};

pub fn expand(input: DeriveInput) -> Result<TokenStream2, syn::Error> {
    let container = ContainerAttrs::parse(&input.attrs)?;

    if container.transparent {
        return expand_transparent(&input);
    }

    match &input.data {
        Data::Struct(DataStruct {
            fields: Fields::Named(named),
            ..
        }) => expand_struct(&input.ident, &named.named, &container),
        Data::Enum(data_enum) => {
            if container.deny_unknown_fields {
                return Err(syn::Error::new_spanned(
                    &input.ident,
                    "deny_unknown_fields is not supported on enums",
                ));
            }
            expand_enum(&input.ident, data_enum)
        }
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
    container: &ContainerAttrs,
) -> Result<TokenStream2, syn::Error> {
    let deny_unknown_fields = container.deny_unknown_fields;
    let mut bindings = Vec::new();
    let mut stream_decls = Vec::new();
    let mut stream_arms = Vec::new();
    let mut stream_finals = Vec::new();
    let mut field_idents = Vec::new();
    let mut known_keys = Vec::new();

    for field in fields {
        let field_info = parse_field_with(field, container.default)?;
        let ident = field_info.ident;
        let key = field_info.key;
        let kind = field_info.kind;

        if matches!(kind, FieldKind::Skip(_)) {
            field_idents.push(ident.clone());
            let ty = match &kind {
                FieldKind::Skip(ty) => ty,
                _ => unreachable!(),
            };
            bindings.push(quote! {
                let #ident: #ty = <#ty as ::std::default::Default>::default();
            });
            stream_decls.push(quote! {
                let #ident: #ty = <#ty as ::std::default::Default>::default();
            });
            continue;
        }

        known_keys.push(key.clone());

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
            FieldKind::Skip(_) => unreachable!(),
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
            FieldKind::Skip(_) => unreachable!(),
        };

        field_idents.push(ident);
        bindings.push(binding);
        stream_decls.push(decl);
        stream_arms.push(arm);
        stream_finals.push(final_bind);
    }

    let finalize = if deny_unknown_fields {
        quote! { h.finalize()?; }
    } else {
        quote! {}
    };

    let stream_unknown_key = if deny_unknown_fields {
        quote! {
            _ => {
                return ::std::result::Result::Err(::umol_edn::EdnError::De(
                    ::umol_edn::DeError::UnknownField {
                        key: __key.into_owned(),
                        path: ::std::vec::Vec::new(),
                    },
                ));
            }
        }
    } else {
        quote! {
            _ => __de.read_skip_value()?,
        }
    };

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
                #finalize
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
                        #stream_unknown_key
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
        let rename = FieldAttrs::parse(&variant.attrs)?.rename;
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
    let mut stream_unit_arms = Vec::new();
    let mut stream_map_arms = Vec::new();

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
                stream_unit_arms.push(quote! {
                    #key => ::std::result::Result::Ok(#name::#ident),
                });
                stream_map_arms.push(quote! {
                    #key => {
                        __de.read_skip_value()?;
                        ::std::result::Result::Ok(#name::#ident)
                    }
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
                stream_map_arms.push(quote! {
                    #key => {
                        let __slice = __de.read_value_slice()?;
                        let __inner = <#inner_ty as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?;
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

                let stream_field_parsers: Vec<_> = fields
                    .unnamed
                    .iter()
                    .map(|f| {
                        let ty = &f.ty;
                        quote! {{
                            let __slice = __de.read_value_slice()?;
                            <#ty as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?
                        }}
                    })
                    .collect();

                stream_map_arms.push(quote! {
                    #key => {
                        __de.consume_byte(b'[')?;
                        #(let __v = #stream_field_parsers;)*
                        // ^ each binding shadows, so collect into tuple below
                        ::std::result::Result::Err(::umol_edn::EdnError::De(
                            ::umol_edn::DeError::Custom("internal".into()),
                        ))
                    }
                });
                // The above doesn't work well for binding N fields. Use
                // read_value_slice on the whole vector and delegate to from_edn_str.
                // Rewrite:
                stream_map_arms.pop();
                stream_map_arms.push(quote! {
                    #key => {
                        let __slice = __de.read_value_slice()?;
                        let __tree = ::umol_edn::read_string(__slice)?;
                        let __seq = match &__tree {
                            ::umol_edn::Edn::Vector(s) => s,
                            other => return ::std::result::Result::Err(
                                ::umol_edn::EdnError::De(::umol_edn::DeError::TypeMismatch {
                                    expected: "vector",
                                    got: other.kind(),
                                    path: ::std::vec::Vec::new(),
                                }),
                            ),
                        };
                        if __seq.len() != #count {
                            return ::std::result::Result::Err(::umol_edn::EdnError::De(
                                ::umol_edn::DeError::Custom(format!(
                                    "expected {} elements for {} variant {}, got {}",
                                    #count, #type_name, #key, __seq.len(),
                                )),
                            ));
                        }
                        ::std::result::Result::Ok(#name::#ident(#(#field_parsers),*))
                    }
                });
            }
            Fields::Named(fields) => {
                let parsed_fields: Vec<_> = fields
                    .named
                    .iter()
                    .map(parse_field)
                    .collect::<Result<_, _>>()?;

                let mut field_bindings = Vec::new();
                let mut stream_decls = Vec::new();
                let mut stream_inner_arms = Vec::new();
                let mut stream_finals = Vec::new();

                for fi in &parsed_fields {
                    let f_ident = &fi.ident;
                    let f_key = &fi.key;
                    match &fi.kind {
                        FieldKind::Skip(ty) => {
                            let b = quote! {
                                let #f_ident: #ty =
                                    <#ty as ::std::default::Default>::default();
                            };
                            field_bindings.push(b.clone());
                            stream_decls.push(b);
                        }
                        FieldKind::Option(inner) => {
                            field_bindings.push(quote! {
                                let #f_ident: ::std::option::Option<#inner> =
                                    __h.optional(#f_key)?;
                            });
                            stream_decls.push(quote! {
                                let mut #f_ident: ::std::option::Option<#inner> =
                                    ::std::option::Option::None;
                            });
                            stream_inner_arms.push(quote! {
                                #f_key => {
                                    let __slice = __de.read_value_slice()?;
                                    #f_ident = ::std::option::Option::Some(
                                        <#inner as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?,
                                    );
                                }
                            });
                        }
                        FieldKind::Defaulted(ty) => {
                            field_bindings.push(quote! {
                                let #f_ident: #ty =
                                    __h.optional::<#ty>(#f_key)?.unwrap_or_default();
                            });
                            stream_decls.push(quote! {
                                let mut #f_ident: #ty =
                                    <#ty as ::std::default::Default>::default();
                            });
                            stream_inner_arms.push(quote! {
                                #f_key => {
                                    let __slice = __de.read_value_slice()?;
                                    #f_ident =
                                        <#ty as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?;
                                }
                            });
                        }
                        FieldKind::Required(ty) => {
                            field_bindings.push(quote! {
                                let #f_ident: #ty = __h.required(#f_key)?;
                            });
                            stream_decls.push(quote! {
                                let mut #f_ident: ::std::option::Option<#ty> =
                                    ::std::option::Option::None;
                            });
                            stream_inner_arms.push(quote! {
                                #f_key => {
                                    let __slice = __de.read_value_slice()?;
                                    #f_ident = ::std::option::Option::Some(
                                        <#ty as ::umol_edn::FromEdn<'de>>::from_edn_str(__slice)?,
                                    );
                                }
                            });
                            stream_finals.push(quote! {
                                let #f_ident = #f_ident.ok_or_else(||
                                    ::umol_edn::EdnError::De(::umol_edn::DeError::MissingField {
                                        key: #f_key.to_string(),
                                        path: ::std::vec::Vec::new(),
                                    })
                                )?;
                            });
                        }
                    }
                }

                let f_idents: Vec<_> = parsed_fields.iter().map(|fi| &fi.ident).collect();

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

                stream_map_arms.push(quote! {
                    #key => {
                        #(#stream_decls)*
                        __de.consume_byte(b'{')?;
                        loop {
                            if __de.try_consume_byte(b'}')? {
                                break;
                            }
                            let __inner_key = __de.read_keyword_name()?;
                            match __inner_key.as_ref() {
                                #(#stream_inner_arms)*
                                _ => __de.read_skip_value()?,
                            }
                        }
                        #(#stream_finals)*
                        ::std::result::Result::Ok(#name::#ident { #(#f_idents),* })
                    }
                });
            }
        }
    }

    let has_unit_variants = !unit_arms.is_empty();
    let has_data_variants = !map_arms.is_empty();

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

    let stream_keyword_branch = if has_unit_variants {
        quote! {
            Some(b':') => {
                let __kw = __de.read_keyword_name()?;
                match __kw.as_ref() {
                    #(#stream_unit_arms)*
                    _ => ::std::result::Result::Err(::umol_edn::EdnError::De(
                        ::umol_edn::DeError::Custom(
                            format!("unknown {} variant: {:?}", #type_name, __kw),
                        ),
                    )),
                }
            }
        }
    } else {
        quote! {}
    };

    let stream_map_branch = if has_data_variants {
        quote! {
            Some(b'{') => {
                __de.consume_byte(b'{')?;
                let __tag = __de.read_keyword_name()?;
                let __result = match __tag.as_ref() {
                    #(#stream_map_arms)*
                    _ => ::std::result::Result::Err(::umol_edn::EdnError::De(
                        ::umol_edn::DeError::Custom(
                            format!("unknown {} variant: {:?}", #type_name, __tag),
                        ),
                    )),
                };
                __de.consume_byte(b'}')?;
                __result
            }
        }
    } else {
        quote! {}
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

            fn from_edn_str(
                input: &'de str,
            ) -> ::std::result::Result<Self, ::umol_edn::EdnError> {
                let mut __de = ::umol_edn::EdnStreamDeserializer::new(input);
                let __result = match __de.peek_byte()? {
                    #stream_keyword_branch
                    #stream_map_branch
                    __other => ::std::result::Result::Err(::umol_edn::EdnError::De(
                        ::umol_edn::DeError::TypeMismatch {
                            expected: #expected,
                            got: match __other {
                                Some(_) => "other",
                                None => "eof",
                            },
                            path: ::std::vec::Vec::new(),
                        },
                    )),
                };
                __de.expect_eof()?;
                __result
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Transparent structs
// ---------------------------------------------------------------------------

fn expand_transparent(input: &DeriveInput) -> Result<TokenStream2, syn::Error> {
    let name = &input.ident;
    let (inner_ty, constructor) = match &input.data {
        Data::Struct(DataStruct { fields: Fields::Named(named), .. }) => {
            if named.named.len() != 1 {
                return Err(syn::Error::new_spanned(
                    name,
                    "transparent requires exactly one field",
                ));
            }
            let field = named.named.first().unwrap();
            let ident = field.ident.as_ref().unwrap();
            (&field.ty, quote! { Self { #ident: __inner } })
        }
        Data::Struct(DataStruct { fields: Fields::Unnamed(unnamed), .. }) => {
            if unnamed.unnamed.len() != 1 {
                return Err(syn::Error::new_spanned(
                    name,
                    "transparent requires exactly one field",
                ));
            }
            (&unnamed.unnamed.first().unwrap().ty, quote! { Self(__inner) })
        }
        _ => {
            return Err(syn::Error::new_spanned(
                name,
                "transparent can only be used on single-field structs",
            ));
        }
    };

    Ok(quote! {
        impl<'de> ::umol_edn::FromEdn<'de> for #name {
            fn from_edn(
                edn: &::umol_edn::Edn<'de>,
            ) -> ::std::result::Result<Self, ::umol_edn::DeError> {
                let __inner = <#inner_ty as ::umol_edn::FromEdn<'de>>::from_edn(edn)?;
                ::std::result::Result::Ok(#constructor)
            }

            fn from_edn_str(
                input: &'de str,
            ) -> ::std::result::Result<Self, ::umol_edn::EdnError> {
                let __inner = <#inner_ty as ::umol_edn::FromEdn<'de>>::from_edn_str(input)?;
                ::std::result::Result::Ok(#constructor)
            }
        }
    })
}

// ---------------------------------------------------------------------------
// Container attributes
// ---------------------------------------------------------------------------

struct ContainerAttrs {
    transparent: bool,
    deny_unknown_fields: bool,
    default: bool,
}

impl ContainerAttrs {
    fn parse(attrs: &[syn::Attribute]) -> Result<Self, syn::Error> {
        let mut transparent = false;
        let mut deny_unknown_fields = false;
        let mut default = false;
        for attr in attrs {
            if !attr.path().is_ident("edn") {
                continue;
            }
            attr.parse_nested_meta(|meta| {
                if meta.path.is_ident("transparent") {
                    transparent = true;
                    Ok(())
                } else if meta.path.is_ident("deny_unknown_fields") {
                    deny_unknown_fields = true;
                    Ok(())
                } else if meta.path.is_ident("default") {
                    default = true;
                    Ok(())
                } else {
                    Err(meta.error("unsupported container #[edn(...)] attribute"))
                }
            })?;
        }
        if transparent && (deny_unknown_fields || default) {
            return Err(syn::Error::new(
                proc_macro2::Span::call_site(),
                "transparent cannot be combined with other container attributes",
            ));
        }
        Ok(ContainerAttrs { transparent, deny_unknown_fields, default })
    }
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
    Skip(Type),
    Option(Type),
    Defaulted(Type),
    Required(Type),
}

fn parse_field(field: &Field) -> Result<FieldInfo, syn::Error> {
    parse_field_with(field, false)
}

fn parse_field_with(field: &Field, container_default: bool) -> Result<FieldInfo, syn::Error> {
    let ident = field
        .ident
        .clone()
        .ok_or_else(|| syn::Error::new_spanned(field, "field must be named"))?;

    let attrs = FieldAttrs::parse(&field.attrs)?;
    let key = attrs.rename.unwrap_or_else(|| ident.to_string().to_kebab_case());

    let kind = if attrs.skip {
        FieldKind::Skip(field.ty.clone())
    } else if let Some(inner) = inner_of(&field.ty, "Option") {
        FieldKind::Option(inner.clone())
    } else if attrs.default || container_default {
        FieldKind::Defaulted(field.ty.clone())
    } else {
        FieldKind::Required(field.ty.clone())
    };

    Ok(FieldInfo { ident, key, kind })
}

struct FieldAttrs {
    rename: Option<String>,
    default: bool,
    skip: bool,
}

impl FieldAttrs {
    fn parse(attrs: &[syn::Attribute]) -> Result<Self, syn::Error> {
        let mut rename = None;
        let mut default = false;
        let mut skip = false;
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
                    default = true;
                    Ok(())
                } else if meta.path.is_ident("skip") {
                    skip = true;
                    Ok(())
                } else if meta.path.is_ident("skip_if") {
                    // skip_if is ser-only, accept silently
                    let _: syn::LitStr = meta.value()?.parse()?;
                    Ok(())
                } else {
                    Err(meta.error("unsupported #[edn(...)] attribute"))
                }
            })?;
        }
        Ok(FieldAttrs { rename, default, skip })
    }
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


