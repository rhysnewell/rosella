use std::collections::BTreeMap;

use super::{ContigMarkers, Place, Tally};
use crate::quality::Quality;

pub const COPY_BINS: usize = 6;

pub const GTDB: &str = "gtdb";
pub const CHECKM: &str = "checkm";

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Reading {
    pub completeness: f64,
    pub contamination: f64,
    pub markers: usize,
    pub groups: usize,
    pub copies: [u32; COPY_BINS],
}

impl Reading {
    fn of(
        markers: usize,
        groups: usize,
        copies: impl Iterator<Item = u32>,
        scored: (f64, f64),
    ) -> Self {
        let mut reading = Self {
            completeness: scored.0,
            contamination: scored.1,
            markers,
            groups,
            ..Default::default()
        };
        for held in copies {
            reading.copies[(held as usize).min(COPY_BINS - 1)] += 1;
        }
        reading
    }

    pub fn quality(&self, set: usize) -> Quality {
        Quality {
            completeness: self.completeness,
            contamination: self.contamination,
            set: set as u16,
        }
    }
}

pub struct SetReading {
    pub set: usize,
    pub gtdb: Reading,
    pub checkm: Option<Reading>,
}

pub struct Duplicate {
    pub panel: &'static str,
    pub marker: String,
    pub contig: usize,
    pub copies: u32,
    pub genes: Vec<Place>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Coding {
    pub genes: u64,
    pub coding_bases: u64,
    pub annotated_bp: u64,
}

// A partial hit marks a marker present but only whole genes count as a second copy.
fn gtdb_copies(tally: &Tally) -> u32 {
    tally.complete.max(u32::from(tally.any > 0))
}

pub(super) struct Doubled {
    pub(super) at: usize,
    pub(super) reads_gtdb: bool,
    pub(super) models: Vec<u16>,
}

impl ContigMarkers {
    pub(super) fn gtdb_on(&self, counts: &[Tally], set: usize) -> Reading {
        let held = counts
            .iter()
            .enumerate()
            .filter(|(marker, _)| self.set.sets.holds(set, *marker))
            .map(|(_, tally)| gtdb_copies(tally))
            .collect::<Vec<_>>();
        if held.is_empty() {
            return Reading::default();
        }
        let total = held.len() as f64;
        let present = held.iter().filter(|copies| **copies > 0).count() as f64;
        let extra = held
            .iter()
            .map(|copies| f64::from(copies.saturating_sub(1)))
            .sum::<f64>();
        Reading::of(
            held.len(),
            held.len(),
            held.iter().copied(),
            (100.0 * present / total, 100.0 * extra / total),
        )
    }

    fn lineage_copies(&self, contigs: &[usize], counts: &[Tally], at: usize) -> Option<Vec<u32>> {
        let searched = contigs
            .iter()
            .map(|contig| self.checkm.get(*contig)?.as_deref())
            .collect::<Option<Vec<_>>>()?;
        let panel = &self.set.checkm;
        let mut copies = vec![0u32; panel.len().max(counts.len())];
        if panel.reads_gtdb(at) {
            for model in panel.models(at) {
                copies[model as usize] = counts.get(model as usize).map_or(0, |tally| tally.any);
            }
            return Some(copies);
        }
        for entry in searched.into_iter().flatten() {
            if usize::from(entry.set) == at {
                copies[entry.model as usize] += u32::from(entry.copies);
            }
        }
        Some(copies)
    }

    pub(super) fn checkm_on(
        &self,
        contigs: &[usize],
        counts: &[Tally],
        set: usize,
    ) -> Option<Reading> {
        let panel = &self.set.checkm;
        let at = panel.lineage(self.set.sets.name(set))?;
        let copies = self.lineage_copies(contigs, counts, at)?;
        let models = panel.models(at);
        Some(Reading::of(
            models.len(),
            panel.group_count(at),
            models.iter().map(|model| copies[*model as usize]),
            panel.score(at, |model| copies[model as usize]),
        ))
    }

    pub fn readings(&self, contigs: &[usize]) -> Option<(usize, Vec<SetReading>)> {
        let (chosen, counts) = self.chosen(contigs)?;
        let sets = (0..self.set.sets.len())
            .map(|set| SetReading {
                set,
                gtdb: self.gtdb_on(&counts, set),
                checkm: self.checkm_on(contigs, &counts, set),
            })
            .collect();
        Some((chosen, sets))
    }

    pub(super) fn doubled(&self, contigs: &[usize]) -> Option<Doubled> {
        let (chosen, counts) = self.chosen(contigs)?;
        let panel = &self.set.checkm;
        let at = panel.lineage(self.set.sets.name(chosen))?;
        let copies = self.lineage_copies(contigs, &counts, at)?;
        Some(Doubled {
            at,
            reads_gtdb: panel.reads_gtdb(at),
            models: panel
                .models(at)
                .into_iter()
                .filter(|model| copies[*model as usize] > 1)
                .collect(),
        })
    }

    pub fn duplicates(&self, contigs: &[usize]) -> Vec<Duplicate> {
        let Some((chosen, counts)) = self.chosen(contigs) else {
            return Vec::new();
        };
        let mut rows = Vec::new();
        for (marker, tally) in counts.iter().enumerate() {
            if self.set.sets.holds(chosen, marker) && gtdb_copies(tally) > 1 {
                rows.extend(self.carrying(GTDB, contigs, marker as u16, true));
            }
        }
        let Some(doubled) = self.doubled(contigs) else {
            return rows;
        };
        for model in doubled.models {
            if doubled.reads_gtdb {
                rows.extend(self.carrying(CHECKM, contigs, model, false));
                continue;
            }
            let mut held = BTreeMap::<usize, u32>::new();
            for contig in contigs {
                let copies = self.checkm[*contig]
                    .iter()
                    .flatten()
                    .filter(|entry| usize::from(entry.set) == doubled.at && entry.model == model);
                for entry in copies {
                    *held.entry(*contig).or_default() += u32::from(entry.copies);
                }
            }
            rows.extend(held.into_iter().map(|(contig, copies)| Duplicate {
                panel: CHECKM,
                marker: self.set.checkm.name(model).to_string(),
                contig,
                copies,
                genes: Vec::new(),
            }));
        }
        rows
    }

    fn carrying<'a>(
        &'a self,
        panel: &'static str,
        contigs: &'a [usize],
        marker: u16,
        whole_only: bool,
    ) -> impl Iterator<Item = Duplicate> + 'a {
        contigs.iter().filter_map(move |contig| {
            let genes = self.per_contig[*contig]
                .iter()
                .filter(|hit| hit.marker == marker && !(whole_only && hit.partial))
                .map(|hit| hit.place)
                .collect::<Vec<_>>();
            (!genes.is_empty()).then(|| Duplicate {
                panel,
                marker: self.set.name(marker).to_string(),
                contig: *contig,
                copies: genes.len() as u32,
                genes,
            })
        })
    }

    pub fn coding(&self, contigs: &[usize]) -> Option<Coding> {
        let mut coding = Coding::default();
        let mut annotated = false;
        for contig in contigs {
            if !self.annotated.get(*contig).copied().unwrap_or_default() {
                continue;
            }
            annotated = true;
            let shape = self.shapes[*contig];
            coding.genes += u64::from(shape.genes);
            coding.coding_bases += shape.coding_bases;
            coding.annotated_bp += self.lengths.get(*contig).copied().unwrap_or_default() as u64;
        }
        annotated.then_some(coding)
    }
}
