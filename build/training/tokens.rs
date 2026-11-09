use std::path::{Path, PathBuf};

use proc_macro2::{Delimiter, Spacing, TokenStream, TokenTree};

const LOGGING: [&str; 5] = ["trace", "debug", "info", "warn", "error"];
const INERT: [&str; 3] = ["doc", "allow", "expect"];
const INCLUDES: [&str; 2] = ["include_str", "include_bytes"];

pub fn is_logging(path: &syn::Path) -> bool {
    let names = path
        .segments
        .iter()
        .map(|segment| segment.ident.to_string())
        .collect::<Vec<_>>();
    match names.as_slice() {
        [name] => LOGGING.contains(&name.as_str()),
        [first, name] => first == "log" && LOGGING.contains(&name.as_str()),
        _ => false,
    }
}

pub fn idents(stream: TokenStream) -> Vec<String> {
    let mut found = Vec::new();
    for tree in stream {
        match tree {
            TokenTree::Ident(ident) => found.push(ident.to_string()),
            TokenTree::Group(group) => found.extend(idents(group.stream())),
            _ => {}
        }
    }
    found
}

pub fn included(stream: TokenStream, file: &Path) -> Vec<PathBuf> {
    let trees = stream.into_iter().collect::<Vec<_>>();
    let mut found = Vec::new();
    for (at, tree) in trees.iter().enumerate() {
        match tree {
            TokenTree::Ident(ident) if INCLUDES.contains(&ident.to_string().as_str()) => {
                let Some(TokenTree::Group(args)) = trees.get(at + 2) else {
                    continue;
                };
                let literal = args.stream().into_iter().find_map(|tree| match tree {
                    TokenTree::Literal(literal) => Some(literal.to_string()),
                    _ => None,
                });
                if let Some(name) = literal.and_then(|l| syn::parse_str::<syn::LitStr>(&l).ok()) {
                    found.push(file.parent().unwrap_or(Path::new(".")).join(name.value()));
                }
            }
            TokenTree::Group(group) => found.extend(included(group.stream(), file)),
            _ => {}
        }
    }
    found
}

// Comments, doc text, lint attributes and log lines leave the trained weights as they were, so
// they leave the key as it was too.
pub fn normalised(stream: TokenStream) -> String {
    let mut out = String::new();
    write_tokens(stream, &mut out);
    out
}

fn write_tokens(stream: TokenStream, out: &mut String) {
    let trees = stream.into_iter().collect::<Vec<_>>();
    let mut at = 0;
    while at < trees.len() {
        if let Some(skip) = skipped(&trees[at..]) {
            at += skip;
            continue;
        }
        match &trees[at] {
            TokenTree::Group(group) => {
                let (open, close) = match group.delimiter() {
                    Delimiter::Parenthesis => ("(", ")"),
                    Delimiter::Brace => ("{", "}"),
                    Delimiter::Bracket => ("[", "]"),
                    Delimiter::None => ("", ""),
                };
                out.push_str(open);
                write_tokens(group.stream(), out);
                out.push_str(close);
            }
            TokenTree::Punct(punct) => {
                out.push(punct.as_char());
                if punct.spacing() == Spacing::Alone {
                    out.push(' ');
                }
            }
            other => {
                out.push_str(&other.to_string());
                out.push(' ');
            }
        }
        at += 1;
    }
}

fn is_punct(tree: Option<&TokenTree>, wanted: char) -> bool {
    matches!(tree, Some(TokenTree::Punct(punct)) if punct.as_char() == wanted)
}

fn skipped(trees: &[TokenTree]) -> Option<usize> {
    if is_punct(trees.first(), '#') {
        let inner = usize::from(is_punct(trees.get(1), '!'));
        return match trees.get(1 + inner) {
            Some(TokenTree::Group(group))
                if group.delimiter() == Delimiter::Bracket
                    && matches!(
                        group.stream().into_iter().next(),
                        Some(TokenTree::Ident(name)) if INERT.contains(&name.to_string().as_str())
                    ) =>
            {
                Some(2 + inner)
            }
            _ => None,
        };
    }
    let path = match trees {
        [TokenTree::Ident(log), first, second, ..]
            if log == "log" && is_punct(Some(first), ':') && is_punct(Some(second), ':') =>
        {
            3
        }
        _ => 0,
    };
    match (trees.get(path), trees.get(path + 1), trees.get(path + 2)) {
        (Some(TokenTree::Ident(name)), bang, Some(TokenTree::Group(_)))
            if LOGGING.contains(&name.to_string().as_str()) && is_punct(bang, '!') =>
        {
            Some(path + 3 + usize::from(is_punct(trees.get(path + 3), ';')))
        }
        _ => None,
    }
}
