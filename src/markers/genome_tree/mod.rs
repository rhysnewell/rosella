use std::path::{Path, PathBuf};

use anyhow::{Result, bail};

use crate::markers::Reading;
use crate::markers::checkm::Panel;

mod fetch;
mod nodes;
mod place;
mod reading;

pub use nodes::{Chosen, Tree, parse_newick};
pub use place::References;

#[derive(Clone, Debug)]
pub struct Lineage {
    pub name: String,
    pub reading: Reading,
}

pub struct Data {
    folder: PathBuf,
    tree: Tree,
    references: References,
    profiles: Vec<place::Profile>,
    phylo: Panel,
    models: String,
    clans: String,
}

impl Data {
    pub fn prepare(directory: &Path) -> Result<Self> {
        let folder = directory.join(fetch::FOLDER);
        if !fetch::complete(&folder) {
            fetch::fetch(directory)?;
        }
        Self::load(&folder)
    }

    fn load(folder: &Path) -> Result<Self> {
        let read = |file: &str| std::fs::read_to_string(folder.join(file));
        let tree = Tree::new(&read(fetch::TREE)?, &read(fetch::NODES)?, fetch::SELECTED)?;
        let references = References::read(&folder.join(fetch::REFERENCES))?;
        let profiles = place::Profile::read(&folder.join(fetch::PHYLO))?;
        let width = profiles.iter().map(|profile| profile.length).sum::<usize>();
        if width != references.width() {
            bail!(
                "the phylogenetic models span {width} columns but the reference alignment {}",
                references.width()
            );
        }
        let models = read(fetch::MODELS)?;
        let clans = read(fetch::CLANS)?;
        let names = profiles.iter().map(|profile| profile.name.as_str());
        let phylo = panel("phylo", names, &models, &clans);
        Ok(Self {
            folder: folder.to_path_buf(),
            tree,
            references,
            profiles,
            phylo,
            models,
            clans,
        })
    }

    fn panel<'a>(&self, name: &str, models: impl IntoIterator<Item = &'a str>) -> Panel {
        panel(name, models, &self.models, &self.clans)
    }
}

// Each bin is scored on its own pruned groups, so the panel needs only the models.
fn panel<'a>(
    name: &str,
    models: impl IntoIterator<Item = &'a str>,
    table: &str,
    clans: &str,
) -> Panel {
    let mut sets = String::from("set\tpanel\tmodel\tgroup\n");
    for model in models {
        sets.push_str(&format!("{name}\tcheckm\t{model}\t0\n"));
    }
    Panel::parse(&sets, table, clans, |_| None)
}
