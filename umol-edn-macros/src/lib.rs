//! Proc macro for constructing `Edn` values from EDN syntax.

use proc_macro2::TokenStream as TokenStream2;
use proc_macro2::{Delimiter, TokenTree};
use quote::quote;

#[proc_macro]
pub fn edn(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    let input2: TokenStream2 = input.into();
    let tokens: Vec<TokenTree> = input2.into_iter().collect();
    match parse_value(&mut tokens.as_slice()) {
        Ok(code) => code.into(),
        Err(msg) => {
            let err = format!("invalid EDN: {msg}");
            quote! { compile_error!(#err) }.into()
        }
    }
}

fn parse_value(tokens: &mut &[TokenTree]) -> Result<TokenStream2, String> {
    skip_discards(tokens)?;
    let tt = tokens.first().ok_or("unexpected end of input")?;
    match tt {
        TokenTree::Literal(lit) => {
            let lit = lit.clone();
            *tokens = &tokens[1..];
            parse_literal(&lit)
        }
        TokenTree::Ident(id) => {
            let id = id.clone();
            *tokens = &tokens[1..];
            parse_ident(&id, tokens)
        }
        TokenTree::Group(g) => {
            let g = g.clone();
            *tokens = &tokens[1..];
            parse_group(&g)
        }
        TokenTree::Punct(_) => parse_punct(tokens),
    }
}

fn skip_discards(tokens: &mut &[TokenTree]) -> Result<(), String> {
    loop {
        if tokens.len() < 2 {
            return Ok(());
        }
        let is_discard = matches!(&tokens[0], TokenTree::Punct(p) if p.as_char() == '#')
            && matches!(&tokens[1], TokenTree::Ident(id) if id.to_string() == "_");
        if !is_discard {
            return Ok(());
        }
        *tokens = &tokens[2..];
        // Consume and discard the next value
        parse_value(tokens)?;
    }
}

fn parse_literal(lit: &proc_macro2::Literal) -> Result<TokenStream2, String> {
    let s = lit.to_string();
    if s.starts_with('"') {
        // String literal — let Rust handle escape parsing
        Ok(quote! { umol_edn::Edn::Str(::std::borrow::Cow::Owned(#lit.to_string())) })
    } else if s.contains('.') || s.contains('e') || s.contains('E') {
        // Float
        let f: f64 = s.parse().map_err(|_| format!("invalid float: {s}"))?;
        Ok(quote! { umol_edn::Edn::Float(#f) })
    } else {
        // Integer — try i64
        let n: i64 = s.parse().map_err(|_| format!("invalid integer: {s}"))?;
        Ok(quote! { umol_edn::Edn::Int(#n) })
    }
}

fn parse_ident(id: &proc_macro2::Ident, tokens: &mut &[TokenTree]) -> Result<TokenStream2, String> {
    let name = id.to_string();
    match name.as_str() {
        "nil" => Ok(quote! { umol_edn::Edn::Nil }),
        "true" => Ok(quote! { umol_edn::Edn::Bool(true) }),
        "false" => Ok(quote! { umol_edn::Edn::Bool(false) }),
        _ => {
            let full = maybe_slashed_name(&name, tokens);
            Ok(quote! { umol_edn::Edn::Symbol(umol_edn::Symbol::new(#full)) })
        }
    }
}

fn parse_group(g: &proc_macro2::Group) -> Result<TokenStream2, String> {
    let inner: Vec<TokenTree> = g.stream().into_iter().collect();
    let mut rest = inner.as_slice();
    match g.delimiter() {
        Delimiter::Bracket => {
            // Vector: [1 2 3]
            let mut elems = Vec::new();
            while !rest.is_empty() {
                elems.push(parse_value(&mut rest)?);
            }
            Ok(quote! { umol_edn::Edn::Vector(vec![#(#elems),*].into()) })
        }
        Delimiter::Parenthesis => {
            // List: (1 2 3)
            let mut elems = Vec::new();
            while !rest.is_empty() {
                elems.push(parse_value(&mut rest)?);
            }
            Ok(quote! { umol_edn::Edn::List(vec![#(#elems),*].into()) })
        }
        Delimiter::Brace => {
            // Map: {:a 1 :b 2}
            parse_map(&mut rest)
        }
        Delimiter::None => Err("unexpected token group".into()),
    }
}

fn parse_map(tokens: &mut &[TokenTree]) -> Result<TokenStream2, String> {
    let mut pairs = Vec::new();
    while !tokens.is_empty() {
        let key = parse_value(tokens)?;
        let val = parse_value(tokens).map_err(|_| "map has odd number of elements")?;
        pairs.push(quote! { map.insert(#key, #val); });
    }
    let len = pairs.len();
    Ok(quote! {{
        let mut map = umol_edn::EdnMap::with_capacity(#len);
        #(#pairs)*
        umol_edn::Edn::Map(map)
    }})
}

fn parse_punct(tokens: &mut &[TokenTree]) -> Result<TokenStream2, String> {
    let p = match &tokens[0] {
        TokenTree::Punct(p) => p.clone(),
        _ => unreachable!(),
    };
    match p.as_char() {
        ':' => {
            // Keyword: :foo or :ns/name
            *tokens = &tokens[1..];
            let name = parse_keyword_name(tokens)?;
            Ok(quote! { umol_edn::Edn::keyword(#name) })
        }
        '#' => {
            *tokens = &tokens[1..];
            parse_hash(tokens)
        }
        '+' | '-' => {
            // Signed number: +1, -3.14
            *tokens = &tokens[1..];
            let sign = p.as_char();
            match tokens.first() {
                Some(TokenTree::Literal(lit)) => {
                    *tokens = &tokens[1..];
                    let s = lit.to_string();
                    if s.contains('.') || s.contains('e') || s.contains('E') {
                        let f: f64 = format!("{sign}{s}")
                            .parse()
                            .map_err(|_| format!("invalid float: {sign}{s}"))?;
                        Ok(quote! { umol_edn::Edn::Float(#f) })
                    } else {
                        let n: i64 = format!("{sign}{s}")
                            .parse()
                            .map_err(|_| format!("invalid integer: {sign}{s}"))?;
                        Ok(quote! { umol_edn::Edn::Int(#n) })
                    }
                }
                Some(TokenTree::Ident(id)) => {
                    // Symbol like -foo or +bar
                    *tokens = &tokens[1..];
                    let name = format!("{sign}{}", id);
                    // Check for ns/name continuation
                    let full = maybe_slashed_name(&name, tokens);
                    Ok(quote! { umol_edn::Edn::Symbol(umol_edn::Symbol::new(#full)) })
                }
                _ => {
                    // Bare symbol / or +
                    let name = sign.to_string();
                    Ok(quote! { umol_edn::Edn::Symbol(umol_edn::Symbol::new(#name)) })
                }
            }
        }
        '/' => {
            // The division symbol
            *tokens = &tokens[1..];
            Ok(quote! { umol_edn::Edn::Symbol(umol_edn::Symbol::new("/")) })
        }
        '\\' => {
            // Character literal: \a, \newline, \u0041
            *tokens = &tokens[1..];
            parse_char_literal(tokens)
        }
        '.' | '*' | '!' | '_' | '?' | '$' | '%' | '&' | '=' | '<' | '>' => {
            // Symbol starting with punctuation
            *tokens = &tokens[1..];
            let start = p.as_char().to_string();
            let name = if let Some(TokenTree::Ident(id)) = tokens.first() {
                let combined = format!("{start}{id}");
                *tokens = &tokens[1..];
                maybe_slashed_name(&combined, tokens)
            } else {
                start
            };
            Ok(quote! { umol_edn::Edn::Symbol(umol_edn::Symbol::new(#name)) })
        }
        _ => Err(format!("unexpected punctuation: '{}'", p.as_char())),
    }
}

fn parse_keyword_name(tokens: &mut &[TokenTree]) -> Result<String, String> {
    match tokens.first() {
        Some(TokenTree::Ident(id)) => {
            *tokens = &tokens[1..];
            let name = id.to_string();
            Ok(maybe_slashed_name(&name, tokens))
        }
        Some(TokenTree::Literal(lit)) => {
            // Digit-start keyword like :0 or :123abc
            *tokens = &tokens[1..];
            let s = lit.to_string();
            Ok(maybe_slashed_name(&s, tokens))
        }
        _ => Err("expected keyword name after ':'".into()),
    }
}

fn maybe_slashed_name(prefix: &str, tokens: &mut &[TokenTree]) -> String {
    // Check for /name continuation
    if tokens.len() >= 2 {
        if let Some(TokenTree::Punct(slash)) = tokens.first() {
            if slash.as_char() == '/' {
                if let Some(TokenTree::Ident(name)) = tokens.get(1) {
                    let full = format!("{prefix}/{name}");
                    *tokens = &tokens[2..];
                    return full;
                }
            }
        }
    }
    prefix.to_string()
}

fn parse_hash(tokens: &mut &[TokenTree]) -> Result<TokenStream2, String> {
    match tokens.first() {
        Some(TokenTree::Group(g)) if g.delimiter() == Delimiter::Brace => {
            // Set: #{ ... }
            *tokens = &tokens[1..];
            let inner: Vec<TokenTree> = g.stream().into_iter().collect();
            let mut rest = inner.as_slice();
            let mut elems = Vec::new();
            while !rest.is_empty() {
                elems.push(parse_value(&mut rest)?);
            }
            Ok(quote! {{
                let mut set = umol_edn::EdnSet::new();
                #(set.insert(#elems);)*
                umol_edn::Edn::Set(set)
            }})
        }
        Some(TokenTree::Punct(p)) if p.as_char() == '#' => {
            Err("## special floats (NaN, Inf, -Inf) are not supported in EDN".into())
        }
        Some(TokenTree::Ident(id)) if id.to_string() == "_" => {
            // #_ discard — already handled in skip_discards, but can appear mid-value
            *tokens = &tokens[1..];
            parse_value(tokens)?; // discard
            parse_value(tokens) // return next
        }
        Some(TokenTree::Ident(id)) => {
            // Tagged literal: #tag value or #ns/tag value
            *tokens = &tokens[1..];
            let tag_name = id.to_string();
            let tag = maybe_slashed_name(&tag_name, tokens);
            let val = parse_value(tokens)?;
            Ok(quote! { umol_edn::Edn::Tagged(#tag.to_string(), Box::new(#val)) })
        }
        _ => Err("unexpected token after #".into()),
    }
}

fn parse_char_literal(tokens: &mut &[TokenTree]) -> Result<TokenStream2, String> {
    match tokens.first() {
        Some(TokenTree::Ident(id)) => {
            *tokens = &tokens[1..];
            let name = id.to_string();
            match name.as_str() {
                "newline" => Ok(quote! { umol_edn::Edn::Char('\n') }),
                "return" => Ok(quote! { umol_edn::Edn::Char('\r') }),
                "space" => Ok(quote! { umol_edn::Edn::Char(' ') }),
                "tab" => Ok(quote! { umol_edn::Edn::Char('\t') }),
                s if s.len() == 1 => {
                    let ch = s.chars().next().unwrap();
                    Ok(quote! { umol_edn::Edn::Char(#ch) })
                }
                s if s.starts_with('u') && s.len() == 5 => {
                    let hex = &s[1..];
                    let cp = u32::from_str_radix(hex, 16)
                        .map_err(|_| format!("invalid unicode escape: \\{s}"))?;
                    let ch = char::from_u32(cp)
                        .ok_or_else(|| format!("invalid code point: \\{s}"))?;
                    Ok(quote! { umol_edn::Edn::Char(#ch) })
                }
                _ => Err(format!("invalid character literal: \\{name}")),
            }
        }
        _ => Err("expected character name after '\\'".into()),
    }
}
