// Procedural macros for umol
#![allow(unused_imports, dead_code)]

use heck::{ToShoutySnakeCase, ToSnakeCase};
use proc_macro::TokenStream;
use proc_macro2::Span;
use quote::{format_ident, quote};
use syn::{
    parse_macro_input, AttributeArgs, GenericArgument, Ident, ItemImpl, Lit, Meta, NestedMeta,
    PathArguments, Type, TypePath,
};

/// Attribute macro that generates model impl code for properties
///
/// This macro must be placed on an `impl Property<M> for P` block. It generates:
/// 1. A constant of the property type for easy access
/// 2. A method on the model type to compute the property
///
/// # Syntax
///
/// ```ignore
/// #[property(method = "atom_count", const_name = "ATOM_COUNT")]
/// impl Property<MyModel> for MyProperty {
///     // ... implementation ...
/// }
/// ```
///
/// or using the alternative syntax:
///
/// ```ignore
/// #[property(method(atom_count), const_name(ATOM_COUNT))]
/// impl Property<MyModel> for MyProperty {
///     // ... implementation ...
/// }
/// ```
///
/// # Parameters
///
/// - `method`: The name of the method to generate on the model type. If not specified,
///   defaults to the snake_case version of the property struct name.
/// - `const_name`: The name of the constant to generate. If not specified, defaults to
///   the SHOUTY_SNAKE_CASE version of the property struct name.
///
/// # Examples
///
/// ```ignore
/// #[property]
/// impl Property<Stoichiometry> for Mass {
///     type Value = f64;
///     type Args = ();
///
///     fn name(&self) -> String {
///         "mass".to_string()
///     }
///
///     fn compute(&self, model: &Stoichiometry, _args: Self::Args) -> Result<Self::Value> {
///         // ... implementation ...
///     }
/// }
/// ```
///
/// This will generate:
///
/// ```ignore
/// impl Stoichiometry {
///     pub const MASS: Mass = Mass;
///
///     pub fn mass(&self) -> Result<f64> {
///         self.compute_property(&Self::MASS, ())
///     }
/// }
/// ```
#[proc_macro_attribute]
pub fn property(attr: TokenStream, item: TokenStream) -> TokenStream {
    // Parse attribute arguments
    let args = parse_macro_input!(attr as AttributeArgs);
    let input = parse_macro_input!(item as ItemImpl);

    // Extract method_name and optional const_name from attributes
    let method_name = extract_method_name(&args).unwrap_or_else(|| {
        // Default to lowercase of property struct name
        let property_struct = extract_property_struct(&input).unwrap();
        format_ident!("{}", property_struct.to_string().to_snake_case())
    });

    let const_name = extract_const_name(&args).unwrap_or_else(|| {
        // Default to uppercased version of property struct name
        let property_struct = extract_property_struct(&input).unwrap();
        format_ident!("{}", property_struct.to_string().to_shouty_snake_case())
    });

    // Extract model and property types from impl
    let model_type = extract_model_type(&input).unwrap();
    let property_struct = extract_property_struct(&input).unwrap();

    // Extract Value and Args types
    let value_type = extract_value_type(&input).unwrap();
    let args_type = extract_args_type(&input).unwrap();

    // Generate the model impl
    let model_impl = if is_unit_type(&args_type) {
        // No args case
        quote! {
            impl #model_type {
                pub const #const_name: #property_struct = #property_struct;

                pub fn #method_name(&self) -> Result<#value_type> {
                    self.compute_property(&Self::#const_name, ())
                }
            }
        }
    } else {
        // With args case
        quote! {
            impl #model_type {
                pub const #const_name: #property_struct = #property_struct;

                pub fn #method_name(&self, args: #args_type) -> Result<#value_type> {
                    self.compute_property(&Self::#const_name, args)
                }
            }
        }
    };

    // Return both the original impl and the new model impl
    let result = quote! {
        #input

        #model_impl
    };

    result.into()
}

// Helper functions
fn extract_method_name(args: &AttributeArgs) -> Option<Ident> {
    for arg in args {
        match arg {
            // Handle #[property(method = "mass")]
            NestedMeta::Meta(Meta::NameValue(nv)) if nv.path.is_ident("method") => {
                if let Lit::Str(lit) = &nv.lit {
                    return Some(Ident::new(&lit.value(), Span::call_site()));
                }
            }
            // Handle #[property(method(mass))]
            NestedMeta::Meta(Meta::List(list)) if list.path.is_ident("method") => {
                if let Some(NestedMeta::Meta(Meta::Path(path))) = list.nested.first() {
                    if let Some(ident) = path.get_ident() {
                        return Some(ident.clone());
                    }
                }
            }
            _ => {}
        }
    }
    None
}

fn extract_const_name(args: &AttributeArgs) -> Option<Ident> {
    for arg in args {
        match arg {
            // Handle #[property(const_name = "MASS")]
            NestedMeta::Meta(Meta::NameValue(nv)) if nv.path.is_ident("const_name") => {
                if let Lit::Str(lit) = &nv.lit {
                    return Some(Ident::new(&lit.value(), Span::call_site()));
                }
            }
            // Handle #[property(const_name(MASS))]
            NestedMeta::Meta(Meta::List(list)) if list.path.is_ident("const_name") => {
                if let Some(NestedMeta::Meta(Meta::Path(path))) = list.nested.first() {
                    if let Some(ident) = path.get_ident() {
                        return Some(ident.clone());
                    }
                }
            }
            _ => {}
        }
    }
    None
}

fn extract_model_type(input: &ItemImpl) -> Option<&Type> {
    // We expect impl Property<ModelType> for PropertyStruct
    if let Some((_, path, _)) = &input.trait_ {
        if path.segments.len() == 1 && path.segments[0].ident == "Property" {
            if let PathArguments::AngleBracketed(args) = &path.segments[0].arguments {
                if let Some(GenericArgument::Type(ty)) = args.args.first() {
                    return Some(ty);
                }
            }
        }
    }
    None
}

fn extract_property_struct(input: &ItemImpl) -> Option<Ident> {
    if let Type::Path(TypePath { path, .. }) = &*input.self_ty {
        if path.segments.len() == 1 {
            return Some(path.segments[0].ident.clone());
        }
    }
    None
}

fn extract_value_type(input: &ItemImpl) -> Option<Type> {
    for item in &input.items {
        if let syn::ImplItem::Type(item_type) = item {
            if item_type.ident == "Value" {
                return Some(item_type.ty.clone());
            }
        }
    }
    None
}

fn extract_args_type(input: &ItemImpl) -> Option<Type> {
    for item in &input.items {
        if let syn::ImplItem::Type(item_type) = item {
            if item_type.ident == "Args" {
                return Some(item_type.ty.clone());
            }
        }
    }
    None
}

fn is_unit_type(ty: &Type) -> bool {
    if let Type::Tuple(tuple) = ty {
        tuple.elems.is_empty()
    } else {
        false
    }
}
