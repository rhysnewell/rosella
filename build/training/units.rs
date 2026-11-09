use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use proc_macro2::TokenStream;
use quote::ToTokens;
use syn::punctuated::Punctuated;
use syn::visit::{self, Visit};

use super::tokens::{idents, included, is_logging, normalised};

#[derive(Default)]
pub struct Refs {
    owner: Option<Vec<String>>,
    pub paths: HashSet<Vec<String>>,
    pub methods: HashSet<String>,
    pub idents: HashSet<String>,
    pub imports: HashMap<String, Vec<String>>,
    pub returns: Vec<Vec<String>>,
}

#[derive(Default)]
pub struct Unit {
    pub file: PathBuf,
    pub module: Vec<String>,
    pub names: Vec<String>,
    pub owner: Option<Vec<String>>,
    pub trait_path: Option<Vec<String>>,
    pub methods: Vec<String>,
    pub import: Option<Vec<String>>,
    pub glob: bool,
    pub shape: bool,
    pub refs: Refs,
    pub text: String,
    pub includes: Vec<PathBuf>,
}

pub struct Place<'a> {
    pub file: &'a Path,
    pub module: &'a [String],
}

impl Place<'_> {
    fn unit(&self, names: Vec<String>, refs: Refs, tokens: TokenStream) -> Unit {
        Unit {
            file: self.file.to_path_buf(),
            module: self.module.to_vec(),
            names,
            refs,
            text: normalised(tokens.clone()),
            includes: included(tokens, self.file),
            ..Unit::default()
        }
    }
}

fn is_test(attrs: &[syn::Attribute]) -> bool {
    attrs.iter().any(|attr| {
        attr.path().is_ident("cfg") && attr.parse_args::<syn::Ident>().is_ok_and(|id| id == "test")
    })
}

fn returned(owner: &Option<Vec<String>>, output: &syn::ReturnType) -> Vec<Vec<String>> {
    let refs = visited(owner.clone(), |refs| refs.visit_return_type(output));
    refs.paths.into_iter().collect()
}

fn visited(owner: Option<Vec<String>>, visit: impl FnOnce(&mut Refs)) -> Refs {
    let mut refs = Refs {
        owner,
        ..Refs::default()
    };
    visit(&mut refs);
    refs
}

pub fn walk(items: &[syn::Item], place: &Place, units: &mut Vec<Unit>) {
    for item in items {
        let unit = |names: Vec<String>| {
            let refs = visited(None, |refs| refs.visit_item(item));
            place.unit(names, refs, item.to_token_stream())
        };
        let named = |ident: &syn::Ident| unit(vec![ident.to_string()]);
        let shaped = |ident: &syn::Ident| Unit {
            shape: true,
            ..named(ident)
        };
        let made = match item {
            syn::Item::Mod(_) => continue,
            syn::Item::Fn(found) if is_test(&found.attrs) => continue,
            syn::Item::Use(found) => {
                for (name, path, glob) in leaves(&found.tree, Vec::new()) {
                    units.push(imported(place, name, path, glob));
                }
                continue;
            }
            syn::Item::Impl(found) => {
                implemented(found, place, units);
                continue;
            }
            syn::Item::Fn(found) => {
                let mut made = named(&found.sig.ident);
                made.refs.returns = returned(&None, &found.sig.output);
                made
            }
            syn::Item::Struct(found) => shaped(&found.ident),
            syn::Item::Enum(found) => shaped(&found.ident),
            syn::Item::Union(found) => shaped(&found.ident),
            syn::Item::Trait(found) => named(&found.ident),
            syn::Item::TraitAlias(found) => named(&found.ident),
            syn::Item::Type(found) => shaped(&found.ident),
            syn::Item::Const(found) => named(&found.ident),
            syn::Item::Static(found) => named(&found.ident),
            syn::Item::ExternCrate(found) => named(&found.ident),
            syn::Item::Macro(found) => match &found.ident {
                Some(ident) => named(ident),
                None => unit(idents(found.mac.tokens.clone())),
            },
            other => unit(idents(other.to_token_stream())),
        };
        units.push(made);
    }
}

fn imported(place: &Place, name: String, path: Vec<String>, glob: bool) -> Unit {
    let text = match glob {
        true => format!("use {}::*", path.join("::")),
        false if path.last() == Some(&name) => format!("use {}", path.join("::")),
        false => format!("use {} as {name}", path.join("::")),
    };
    let mut refs = Refs::default();
    if !glob {
        refs.paths.insert(path.clone());
    }
    Unit {
        text,
        import: Some(path),
        glob,
        refs,
        ..place.unit(vec![name], Refs::default(), TokenStream::new())
    }
}

pub fn leaves(tree: &syn::UseTree, prefix: Vec<String>) -> Vec<(String, Vec<String>, bool)> {
    let with = |name: &syn::Ident| {
        let mut path = prefix.clone();
        path.push(name.to_string());
        path
    };
    match tree {
        syn::UseTree::Path(path) => leaves(&path.tree, with(&path.ident)),
        syn::UseTree::Name(_) if prefix == ["log"] => Vec::new(),
        syn::UseTree::Name(name) if name.ident == "self" => {
            let last = prefix.last().cloned().unwrap_or_default();
            vec![(last, prefix, false)]
        }
        syn::UseTree::Name(name) => vec![(name.ident.to_string(), with(&name.ident), false)],
        syn::UseTree::Rename(rename) => {
            vec![(rename.rename.to_string(), with(&rename.ident), false)]
        }
        syn::UseTree::Glob(_) => vec![("*".to_string(), prefix, true)],
        syn::UseTree::Group(group) => group
            .items
            .iter()
            .flat_map(|tree| leaves(tree, prefix.clone()))
            .collect(),
    }
}

fn segments(path: &syn::Path) -> Vec<String> {
    path.segments.iter().map(|s| s.ident.to_string()).collect()
}

fn type_path(ty: &syn::Type) -> Vec<String> {
    match ty {
        syn::Type::Path(path) => segments(&path.path),
        syn::Type::Reference(inner) => type_path(&inner.elem),
        syn::Type::Slice(inner) => type_path(&inner.elem),
        syn::Type::Array(inner) => type_path(&inner.elem),
        _ => Vec::new(),
    }
}

fn item_names(item: &syn::ImplItem) -> Vec<String> {
    match item {
        syn::ImplItem::Fn(found) => vec![found.sig.ident.to_string()],
        syn::ImplItem::Const(found) => vec![found.ident.to_string()],
        syn::ImplItem::Type(found) => vec![found.ident.to_string()],
        syn::ImplItem::Macro(found) => idents(found.mac.tokens.clone()),
        _ => Vec::new(),
    }
}

// A trait impl runs without its name being written, as Drop and Iterator do, so it is kept whole
// with its type. An inherent method is reached on its own.
fn implemented(found: &syn::ItemImpl, place: &Place, units: &mut Vec<Unit>) {
    let owner = type_path(&found.self_ty);
    if let Some((_, path, _)) = &found.trait_ {
        let mut refs = visited(Some(owner.clone()), |refs| refs.visit_item_impl(found));
        refs.returns = found
            .items
            .iter()
            .filter_map(|item| match item {
                syn::ImplItem::Fn(method) => {
                    Some(returned(&Some(owner.clone()), &method.sig.output))
                }
                _ => None,
            })
            .flatten()
            .collect();
        units.push(Unit {
            owner: Some(owner.clone()),
            trait_path: Some(segments(path)),
            methods: found.items.iter().flat_map(item_names).collect(),
            ..place.unit(Vec::new(), refs, found.to_token_stream())
        });
        return;
    }
    for item in &found.items {
        let mut refs = visited(Some(owner.clone()), |refs| {
            refs.visit_generics(&found.generics);
            refs.visit_type(&found.self_ty);
            refs.visit_impl_item(item);
        });
        if let syn::ImplItem::Fn(method) = item {
            refs.returns = returned(&Some(owner.clone()), &method.sig.output);
        }
        let mut tokens = found.generics.to_token_stream();
        tokens.extend(found.self_ty.to_token_stream());
        tokens.extend(item.to_token_stream());
        units.push(Unit {
            owner: Some(owner.clone()),
            methods: item_names(item),
            ..place.unit(Vec::new(), refs, tokens)
        });
    }
}

impl<'ast> Visit<'ast> for Refs {
    fn visit_path(&mut self, path: &'ast syn::Path) {
        let names = path
            .segments
            .iter()
            .flat_map(|segment| match (segment.ident == "Self", &self.owner) {
                (true, Some(owner)) => owner.clone(),
                _ => vec![segment.ident.to_string()],
            })
            .collect();
        self.paths.insert(names);
        visit::visit_path(self, path);
    }

    fn visit_expr_method_call(&mut self, call: &'ast syn::ExprMethodCall) {
        self.methods.insert(call.method.to_string());
        visit::visit_expr_method_call(self, call);
    }

    fn visit_item_use(&mut self, found: &'ast syn::ItemUse) {
        for (name, path, _) in leaves(&found.tree, Vec::new()) {
            self.imports.insert(name, path);
        }
    }

    fn visit_macro(&mut self, found: &'ast syn::Macro) {
        if is_logging(&found.path) {
            return;
        }
        self.visit_path(&found.path);
        match found.parse_body_with(Punctuated::<syn::Expr, syn::Token![,]>::parse_terminated) {
            Ok(arguments) => arguments
                .iter()
                .for_each(|argument| self.visit_expr(argument)),
            Err(_) => self.idents.extend(idents(found.tokens.clone())),
        }
    }
}
