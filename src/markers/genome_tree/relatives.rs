use std::collections::{BTreeSet, HashMap};

use crate::markers::checkm::bare;

pub struct Relatives {
    clan_of: HashMap<String, String>,
    members: HashMap<String, Vec<String>>,
}

impl Relatives {
    pub fn parse(models: &str, clans: &str) -> Self {
        let clan_of = clans
            .lines()
            .skip(1)
            .filter_map(|line| {
                let mut fields = line.split('\t');
                let pfam = fields.next()?;
                let clan = fields.next().filter(|clan| !clan.is_empty())?;
                Some((pfam.to_string(), clan.to_string()))
            })
            .collect::<HashMap<_, _>>();
        let mut members = HashMap::<String, Vec<String>>::new();
        for model in models
            .lines()
            .skip(1)
            .filter_map(|line| line.split('\t').next())
        {
            if let Some(clan) = clan_of.get(bare(model)) {
                members
                    .entry(clan.clone())
                    .or_default()
                    .push(model.to_string());
            }
        }
        Self { clan_of, members }
    }

    // CheckM1 searches far past the bin's sets, so a stronger relative from the same clan can
    // claim a gene. Only a relative can, so the clans stand in for the rest.
    pub fn widened<'a>(&self, chosen: impl IntoIterator<Item = &'a String>) -> BTreeSet<String> {
        let mut held = BTreeSet::new();
        let mut clans = BTreeSet::new();
        for model in chosen {
            clans.extend(self.clan_of.get(bare(model)));
            held.insert(model.clone());
        }
        for clan in clans {
            held.extend(self.members.get(clan).into_iter().flatten().cloned());
        }
        held
    }
}
