use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

use super::units::Unit;

const IMPLICIT: [&str; 26] = [
    "Drop",
    "Deref",
    "DerefMut",
    "Iterator",
    "IntoIterator",
    "Add",
    "Sub",
    "Mul",
    "Div",
    "Neg",
    "Index",
    "IndexMut",
    "AddAssign",
    "SubAssign",
    "MulAssign",
    "DivAssign",
    "PartialEq",
    "PartialOrd",
    "Ord",
    "Hash",
    "Hasher",
    "From",
    "Default",
    "Sum",
    "FromIterator",
    "Extend",
];

struct Resolver<'a> {
    units: &'a [Unit],
    modules: HashMap<&'a [String], &'a Path>,
    defined: HashMap<(&'a Path, &'a str), Vec<usize>>,
    imports: HashMap<(&'a Path, &'a str), Vec<usize>>,
    globs: HashMap<&'a Path, Vec<usize>>,
    macros: HashMap<&'a str, Vec<usize>>,
}

impl<'a> Resolver<'a> {
    fn new(units: &'a [Unit], files: &'a [(PathBuf, Vec<String>)]) -> Self {
        let mut resolver = Self {
            units,
            modules: files
                .iter()
                .map(|(file, module)| (module.as_slice(), file.as_path()))
                .collect(),
            defined: HashMap::new(),
            imports: HashMap::new(),
            globs: HashMap::new(),
            macros: HashMap::new(),
        };
        for (at, unit) in units.iter().enumerate() {
            let file = unit.file.as_path();
            for name in &unit.names {
                let slot = match (&unit.import, unit.glob) {
                    (Some(_), true) => resolver.globs.entry(file).or_default(),
                    (Some(_), false) => resolver.imports.entry((file, name)).or_default(),
                    (None, _) => resolver.defined.entry((file, name)).or_default(),
                };
                slot.push(at);
                if unit.import.is_none() && unit.text.starts_with("macro_rules") {
                    resolver.macros.entry(name).or_default().push(at);
                }
            }
        }
        resolver
    }

    fn absolute(&self, module: &[String], path: &[String]) -> Option<Vec<String>> {
        let (first, rest) = path.split_first()?;
        match first.as_str() {
            "crate" => Some(rest.to_vec()),
            "self" => Some([module, rest].concat()),
            "super" => {
                let ups = path
                    .iter()
                    .take_while(|segment| *segment == "super")
                    .count();
                let kept = module.len().saturating_sub(ups);
                Some([&module[..kept], &path[ups..]].concat())
            }
            _ => {
                let child = [module, &path[..1]].concat();
                self.modules
                    .contains_key(child.as_slice())
                    .then(|| [module, path].concat())
            }
        }
    }

    fn in_file(&self, file: &Path, name: &str, depth: usize) -> Vec<usize> {
        let mut found = Vec::new();
        for map in [&self.defined, &self.imports] {
            found.extend(map.get(&(file, name)).into_iter().flatten());
        }
        if !found.is_empty() || depth > 8 {
            return found;
        }
        for glob in self.globs.get(file).into_iter().flatten() {
            let unit = &self.units[*glob];
            let module = unit.import.as_deref().unwrap_or_default();
            let Some(path) = self.absolute(&unit.module, module) else {
                continue;
            };
            let through = self.lookup(&[path, vec![name.to_string()]].concat(), depth + 1);
            if !through.is_empty() {
                found.push(*glob);
                found.extend(through);
            }
        }
        found
    }

    fn lookup(&self, path: &[String], depth: usize) -> Vec<usize> {
        let found = (0..=path.len())
            .rev()
            .find_map(|at| Some((*self.modules.get(&path[..at])?, &path[at..])));
        match found {
            Some((file, [name, ..])) => self.in_file(file, name, depth),
            _ => Vec::new(),
        }
    }

    fn resolve(&self, unit: &Unit, path: &[String]) -> Vec<usize> {
        let Some((first, rest)) = path.split_first() else {
            return Vec::new();
        };
        let file = unit.file.as_path();
        if let Some(local) = unit.refs.imports.get(first) {
            return self
                .absolute(&unit.module, local)
                .map(|local| self.lookup(&[local.as_slice(), rest].concat(), 0))
                .unwrap_or_default();
        }
        if rest.is_empty() {
            let found = self.in_file(file, first, 0);
            return match found.is_empty() {
                true => self.macros.get(first.as_str()).cloned().unwrap_or_default(),
                false => found,
            };
        }
        if let Some(absolute) = self.absolute(&unit.module, path) {
            return self.lookup(&absolute, 0);
        }
        let mut found = Vec::new();
        for import in self
            .imports
            .get(&(file, first.as_str()))
            .into_iter()
            .flatten()
        {
            let leaf = &self.units[*import];
            found.push(*import);
            let through = leaf
                .import
                .as_deref()
                .and_then(|path| self.absolute(&leaf.module, path))
                .map(|path| self.lookup(&[path.as_slice(), rest].concat(), 0));
            found.extend(through.into_iter().flatten());
        }
        match found.is_empty() {
            true => self.in_file(file, first, 0),
            false => found,
        }
    }

    fn defs(&self, unit: &Unit, path: &[String], depth: usize) -> Vec<usize> {
        let mut out = Vec::new();
        for at in self.resolve(unit, path) {
            let found = &self.units[at];
            match (&found.import, found.glob) {
                (None, _) => out.push(at),
                (Some(import), false) if depth < 8 => {
                    out.extend(self.defs(found, import, depth + 1))
                }
                _ => {}
            }
        }
        out
    }
}

#[derive(Default)]
struct Seen<'a> {
    named: HashMap<&'a Path, HashSet<usize>>,
    coded: HashSet<usize>,
    shaped: HashSet<usize>,
    called: HashMap<&'a Path, HashSet<String>>,
    qualified: HashSet<(usize, String)>,
}

struct Owned {
    at: usize,
    owners: Vec<usize>,
    traits: Vec<usize>,
    implicit: bool,
}

// A method counts once one reached file names its type and calls it. A field's type alone counts
// only for traits that run unasked, and another crate's trait may run from inside that crate.
impl<'a> Seen<'a> {
    fn holds(&self, owned: &Owned, unit: &Unit) -> bool {
        let any = |found: &HashSet<usize>| {
            owned.owners.is_empty() || owned.owners.iter().any(|at| found.contains(at))
        };
        let in_code = any(&self.coded);
        let asked = |item: &String| {
            owned
                .owners
                .iter()
                .chain(&owned.traits)
                .any(|at| self.qualified.contains(&(*at, item.clone())))
        };
        let called = |item: &String| {
            self.called
                .iter()
                .any(|(file, names)| names.contains(item) && self.named.get(file).is_some_and(any))
        };
        match &unit.trait_path {
            Some(_) if owned.implicit => {
                in_code || owned.owners.iter().any(|at| self.shaped.contains(at))
            }
            Some(_) if owned.traits.is_empty() => in_code,
            Some(_) => {
                in_code
                    && (owned.traits.iter().any(|at| self.coded.contains(at))
                        || unit.methods.iter().any(|item| called(item) || asked(item)))
            }
            None => unit.methods.iter().any(|item| asked(item) || called(item)),
        }
    }

    fn name(&mut self, file: &'a Path, found: &[usize], fields: &[Vec<usize>]) -> bool {
        let named = self.named.entry(file).or_default();
        let mut pending = found.to_vec();
        let mut grew = false;
        while let Some(at) = pending.pop() {
            if named.insert(at) {
                grew = true;
                pending.extend(&fields[at]);
            }
        }
        grew
    }
}

fn picked(units: &[Unit], resolver: &Resolver, (file, name): &(&str, &str)) -> Vec<usize> {
    let found = match *name {
        "*" => units
            .iter()
            .enumerate()
            .filter(|(_, unit)| unit.file.starts_with(file))
            .map(|(at, _)| at)
            .collect(),
        _ => resolver
            .defined
            .get(&(Path::new(file), *name))
            .cloned()
            .unwrap_or_default(),
    };
    assert!(!found.is_empty(), "training item {file}:{name} is gone");
    found
}

pub fn reach(
    units: &[Unit],
    files: &[(PathBuf, Vec<String>)],
    entries: &[(&str, &str)],
    outside: &[(&str, &str)],
    whole: Vec<usize>,
) -> BTreeSet<usize> {
    let resolver = Resolver::new(units, files);
    let mut pending = whole;
    pending.extend(
        entries
            .iter()
            .flat_map(|item| picked(units, &resolver, item)),
    );
    let stops = outside
        .iter()
        .flat_map(|item| picked(units, &resolver, item))
        .collect::<BTreeSet<_>>();
    let owned = units
        .iter()
        .enumerate()
        .filter_map(|(at, unit)| {
            let owner = unit.owner.as_deref()?;
            let traits = unit.trait_path.as_deref().unwrap_or_default();
            Some(Owned {
                at,
                owners: resolver.defs(unit, owner, 0),
                traits: resolver.defs(unit, traits, 0),
                implicit: traits
                    .last()
                    .is_some_and(|name| IMPLICIT.contains(&name.as_str())),
            })
        })
        .collect::<Vec<_>>();
    let resolved = |paths: &[Vec<String>], unit: &Unit| {
        paths
            .iter()
            .flat_map(|path| resolver.defs(unit, path, 0))
            .collect::<Vec<_>>()
    };
    let returned = units
        .iter()
        .map(|unit| resolved(&unit.refs.returns, unit))
        .collect::<Vec<_>>();
    let fields = units
        .iter()
        .map(|unit| match unit.shape {
            true => resolved(&unit.refs.paths.iter().cloned().collect::<Vec<_>>(), unit),
            false => Vec::new(),
        })
        .collect::<Vec<_>>();
    let mut reached = stops.clone();
    let mut seen = Seen::default();
    loop {
        while let Some(at) = pending.pop() {
            if !reached.insert(at) {
                continue;
            }
            let unit = &units[at];
            let file = unit.file.as_path();
            let words = unit.refs.idents.iter().map(|name| vec![name.clone()]);
            for path in unit.refs.paths.iter().cloned().chain(words) {
                pending.extend(resolver.resolve(unit, &path));
                let mut defs = resolver.defs(unit, &path, 0);
                if let [_, .., item] = path.as_slice() {
                    seen.qualified
                        .extend(defs.iter().map(|at| (*at, item.clone())));
                }
                if unit.shape {
                    seen.shaped.extend(defs);
                    continue;
                }
                defs.extend(defs.clone().iter().flat_map(|at| &returned[*at]));
                seen.name(file, &defs, &fields);
                seen.coded.extend(defs);
            }
            let called = seen.called.entry(file).or_default();
            called.extend(unit.refs.methods.iter().cloned());
            called.extend(unit.refs.idents.iter().cloned());
        }
        let mut grew = false;
        for found in owned.iter().filter(|found| reached.contains(&found.at)) {
            let unit = &units[found.at];
            if stops.contains(&found.at) || returned[found.at].is_empty() {
                continue;
            }
            let callers = seen
                .called
                .iter()
                .filter(|(_, names)| unit.methods.iter().any(|item| names.contains(item)))
                .map(|(file, _)| *file)
                .collect::<Vec<_>>();
            for file in callers {
                grew |= seen.name(file, &returned[found.at], &fields);
                seen.coded.extend(&returned[found.at]);
            }
        }
        pending = owned
            .iter()
            .filter(|found| !reached.contains(&found.at) && seen.holds(found, &units[found.at]))
            .map(|found| found.at)
            .collect();
        if pending.is_empty() && !grew {
            break;
        }
    }
    let files = reached
        .iter()
        .filter(|at| !stops.contains(at))
        .map(|at| units[*at].file.as_path())
        .collect::<HashSet<_>>();
    let globs = resolver
        .globs
        .iter()
        .filter(|(file, _)| files.contains(*file))
        .flat_map(|(_, found)| found.iter().copied())
        .collect::<Vec<_>>();
    reached.extend(globs);
    &reached - &stops
}
