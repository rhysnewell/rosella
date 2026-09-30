use std::cmp::Reverse;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::quality::Scorer;

struct Event {
    stage: String,
    core: Vec<usize>,
    cut: Vec<usize>,
}

#[derive(Default)]
pub struct CutLog {
    events: Vec<Event>,
}

pub fn owners<'a, B, C>(bins: B) -> HashMap<usize, usize>
where
    B: IntoIterator<Item = (usize, C)>,
    C: IntoIterator<Item = &'a usize>,
{
    bins.into_iter()
        .flat_map(|(label, members)| members.into_iter().map(move |contig| (*contig, label)))
        .collect()
}

fn heir(
    members: &[usize],
    owner: &HashMap<usize, usize>,
    length: impl Fn(usize) -> usize,
) -> Option<usize> {
    let mut held = HashMap::<usize, usize>::new();
    for contig in members {
        if let Some(bin) = owner.get(contig) {
            *held.entry(*bin).or_default() += length(*contig);
        }
    }
    held.into_iter()
        .max_by_key(|(bin, bp)| (*bp, Reverse(*bin)))
        .map(|(bin, _)| bin)
}

impl CutLog {
    pub fn events(&self) -> impl Iterator<Item = (&str, &[usize], &[usize])> {
        self.events.iter().map(|event| {
            (
                event.stage.as_str(),
                event.core.as_slice(),
                event.cut.as_slice(),
            )
        })
    }

    pub fn diff<'a, B, C>(
        &mut self,
        stage: &str,
        before: B,
        after: &HashMap<usize, usize>,
        length: impl Fn(usize) -> usize + Copy,
    ) where
        B: IntoIterator<Item = C>,
        C: IntoIterator<Item = &'a usize>,
    {
        for members in before {
            let members = members.into_iter().copied().collect::<Vec<_>>();
            self.cut(stage, &members, after, length);
        }
    }

    pub fn cut(
        &mut self,
        stage: &str,
        members: &[usize],
        after: &HashMap<usize, usize>,
        length: impl Fn(usize) -> usize,
    ) {
        let heir = heir(members, after, length);
        let (core, cut): (Vec<usize>, Vec<usize>) = members
            .iter()
            .partition(|contig| heir.is_some() && after.get(contig) == heir.as_ref());
        if !cut.is_empty() {
            self.events.push(Event {
                stage: stage.to_string(),
                core,
                cut,
            });
        }
    }

    // A give-back could only hand a cut to the final bin holding most of what its old bin kept.
    pub fn write(
        &self,
        path: &Path,
        bins: &HashMap<usize, Vec<usize>>,
        quality: &dyn Scorer,
        worth: f64,
        lengths: &[usize],
        names: &[String],
    ) -> Result<()> {
        let owner = owners(bins.iter().map(|(label, members)| (*label, members)));
        let mut base = HashMap::new();
        let mut out = BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            out,
            "stage\tcontig\tlength\tmarkers\tcore_bp\tfinal\tended\tback\t\
             completeness\tcontamination\tcompleteness_with\tcontamination_with\tgain"
        )?;
        for (stage, core, cut) in self.events() {
            let target = heir(core, &owner, |contig| lengths[contig]);
            let core_bp = core.iter().map(|contig| lengths[*contig]).sum::<usize>();
            for contig in cut {
                let ended = owner.get(contig).copied();
                let label = |bin: Option<usize>| bin.map_or("-".to_string(), |bin| bin.to_string());
                let markers = quality.features(&[*contig]).len();
                write!(
                    out,
                    "{stage}\t{}\t{}\t{markers}\t{core_bp}\t{}\t{}",
                    names[*contig],
                    lengths[*contig],
                    label(target),
                    label(ended)
                )?;
                let Some(bin) = target.filter(|bin| ended != Some(*bin)) else {
                    writeln!(out, "\t{}\t-\t-\t-\t-\t-", u8::from(target.is_some()))?;
                    continue;
                };
                let members = &bins[&bin];
                let before = *base.entry(bin).or_insert_with(|| quality.score(members));
                let mut offered = members.clone();
                offered.push(*contig);
                let with = quality.score(&offered);
                writeln!(
                    out,
                    "\t0\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.4}",
                    before.completeness,
                    before.contamination,
                    with.completeness,
                    with.contamination,
                    with.score(worth) - before.score(worth)
                )?;
            }
        }
        out.flush()?;
        Ok(())
    }
}
