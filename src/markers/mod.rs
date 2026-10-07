use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::quality::Quality;

mod annotation;
mod annotator;
pub mod cache;
pub mod checkm;
pub mod eukaryote;
pub mod fragments;
mod hit;
pub mod hmm_table;
mod multicopy;
mod proteins;
mod readings;
pub mod replicon;
pub mod sets;
mod strain;
mod table;
mod walk;

pub use annotation::MarkerAnnotation;
use annotation::inflate;
pub use annotator::Annotator;
pub use hit::{Hit, Place};
pub use readings::{CHECKM, COPY_BINS, Coding, Duplicate, GTDB, Reading, SetReading};
pub use strain::{heterogeneity, identity, rebuild};
pub use table::MarkerSet;
pub use walk::Shed;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");

// A marker gene cut by a contig end is still that marker gene, but two halves of one gene on
// two contigs are not two copies, so presence and duplication read different columns.
#[derive(Clone, Copy, Debug)]
pub struct MarkerRules {
    pub fragment_span: f64,
}

impl Default for MarkerRules {
    fn default() -> Self {
        Self {
            fragment_span: fragments::DEFAULT_SPAN,
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Tally {
    complete: u32,
    any: u32,
}

// The search returns its hits in hash order, and the cache and the report are written from
// these, so two annotations of one assembly would not diff against each other.
fn in_marker_order(per_contig: &mut [Vec<Hit>]) {
    for hits in per_contig {
        hits.sort_unstable_by_key(|hit| {
            (
                hit.marker,
                hit.partial,
                hit.place.gene_begin,
                hit.place.reverse,
            )
        });
    }
}

pub struct ContigMarkers {
    per_contig: Vec<Vec<Hit>>,
    checkm: Vec<Option<Vec<checkm::Copies>>>,
    proteins: Vec<proteins::Proteins>,
    shapes: Vec<replicon::Shape>,
    annotated: Vec<bool>,
    lengths: Vec<usize>,
    set: MarkerSet,
    multicopy: std::sync::OnceLock<Vec<bool>>,
}

impl crate::quality::Scorer for ContigMarkers {
    fn features(&self, contigs: &[usize]) -> std::collections::HashSet<u32> {
        self.counts(contigs)
            .iter()
            .enumerate()
            .filter(|(_, tally)| tally.any > 0)
            .map(|(marker, _)| marker as u32)
            .collect()
    }

    // Read against whichever lineage the bin's pattern of absences fits, since a reduced
    // genome is missing markers a whole one of another lineage would carry.
    fn score(&self, contigs: &[usize]) -> Quality {
        let Some((chosen, counts)) = self.chosen(contigs) else {
            return Quality::default();
        };
        match self.gtdb_on(&counts, chosen) {
            reading if reading.markers == 0 => Quality::default(),
            reading => reading.quality(chosen),
        }
    }
}

fn observed(counts: &[Tally]) -> Vec<u16> {
    counts
        .iter()
        .enumerate()
        .filter(|(_, tally)| tally.any > 0)
        .map(|(marker, _)| marker as u16)
        .collect()
}

impl ContigMarkers {
    pub fn new(mut per_contig: Vec<Vec<Hit>>, set: MarkerSet) -> Self {
        in_marker_order(&mut per_contig);
        Self {
            multicopy: std::sync::OnceLock::new(),
            per_contig,
            checkm: Vec::new(),
            proteins: Vec::new(),
            shapes: Vec::new(),
            annotated: Vec::new(),
            lengths: Vec::new(),
            set,
        }
    }

    pub fn set_name(&self, set: u16) -> &str {
        self.set.sets.name(set as usize)
    }

    pub fn report<'a>(
        &self,
        placed: impl IntoIterator<Item = (&'a str, usize)>,
        names: &[String],
        path: &Path,
    ) -> Result<()> {
        let mut sink = BufWriter::new(std::fs::File::create(path)?);
        writeln!(
            sink,
            "bin\tcontig\tmodel\tpartial\tmodel_from\tmodel_to\tmodel_length\tprotein_from\t\
             protein_to\tscore\tgene_begin\tgene_end\tstrand\tcut_left\tcut_right"
        )?;
        for (bin, contig) in placed {
            for hit in &self.per_contig[contig] {
                writeln!(
                    sink,
                    "{bin}\t{}\t{}\t{}\t{}",
                    names[contig],
                    self.set.name(hit.marker),
                    u8::from(hit.partial),
                    hit.fields('\t')
                )?;
            }
        }
        sink.flush()?;
        Ok(())
    }

    pub fn with_lengths(mut self, lengths: Vec<usize>) -> Self {
        self.lengths = lengths;
        self
    }

    pub fn with_shapes(mut self, shapes: Vec<Option<replicon::Shape>>) -> Self {
        self.annotated = shapes.iter().map(Option::is_some).collect();
        self.shapes = shapes.into_iter().map(Option::unwrap_or_default).collect();
        self
    }

    pub fn with_checkm(mut self, checkm: Vec<Option<Vec<checkm::Copies>>>) -> Self {
        self.checkm = checkm;
        self
    }

    fn with_proteins(mut self, kept: Option<proteins::Proteins>) -> Self {
        self.proteins.extend(kept);
        self
    }

    pub fn departures<'a>(
        &self,
        bins: impl IntoIterator<Item = &'a [usize]>,
        composition: ndarray::ArrayView2<f64>,
        depths: ndarray::ArrayView2<f64>,
    ) -> replicon::Departures {
        if self.shapes.len() != self.lengths.len() {
            return replicon::Departures::default();
        }
        let carries = self
            .per_contig
            .iter()
            .map(|hits| !hits.is_empty())
            .collect::<Vec<_>>();
        replicon::departures(
            &self.shapes,
            &carries,
            &self.lengths,
            bins,
            composition,
            depths,
        )
    }

    fn chosen(&self, contigs: &[usize]) -> Option<(usize, Vec<Tally>)> {
        let counts = self.counts(contigs);
        let chosen = self
            .set
            .sets
            .choose(&observed(&counts), self.bin_bp(contigs))?;
        Some((chosen, counts))
    }

    pub fn points(&self, contigs: &[usize], weight: f64) -> Vec<(usize, f64)> {
        let Some((chosen, counts)) = self.chosen(contigs) else {
            return Vec::new();
        };
        counts
            .iter()
            .enumerate()
            .filter(|(marker, _)| self.set.sets.holds(chosen, *marker))
            .map(|(marker, tally)| {
                let extra = f64::from(tally.complete.saturating_sub(1));
                let present = f64::from(u8::from(tally.any >= 1));
                (marker, 100.0 * (present - weight * extra))
            })
            .collect()
    }

    pub fn hit_count(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .map(|contig| self.per_contig[*contig].len())
            .sum()
    }

    fn bin_bp(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .filter_map(|contig| self.lengths.get(*contig))
            .sum()
    }

    // Whether the contig carries a whole marker copy the rest of the bin does not, which is
    // the difference between a home and a second copy of what is already there.
    pub fn completes(&self, contigs: &[usize], contig: usize) -> bool {
        self.repeats(contigs, contig) == Some(false)
    }

    // None when the contig holds no whole marker of the bin's set, so it is no evidence either way.
    pub fn repeats(&self, contigs: &[usize], contig: usize) -> Option<bool> {
        self.repeated(contigs, contig, true, |_, _| true)
    }

    pub fn repeats_any(&self, contigs: &[usize], contig: usize) -> Option<bool> {
        self.repeated(contigs, contig, false, |_, _| true)
    }

    pub fn repeats_in_place(&self, contigs: &[usize], contig: usize) -> Option<bool> {
        self.repeated(contigs, contig, false, Hit::same_part)
    }

    fn repeated(
        &self,
        contigs: &[usize],
        contig: usize,
        whole_only: bool,
        same: impl Fn(&Hit, &Hit) -> bool,
    ) -> Option<bool> {
        let (chosen, _) = self.chosen(contigs)?;
        let counted = |hit: &&Hit| {
            !(whole_only && hit.partial)
                && self.set.sets.holds(chosen, hit.marker as usize)
                && !self.is_multicopy(hit.marker)
        };
        let held = self.per_contig[contig]
            .iter()
            .filter(counted)
            .collect::<Vec<_>>();
        if held.is_empty() {
            return None;
        }
        let rest = contigs
            .iter()
            .filter(|member| **member != contig)
            .flat_map(|member| self.per_contig[*member].iter().filter(counted))
            .filter(|other| held.iter().any(|hit| hit.marker == other.marker))
            .collect::<Vec<_>>();
        Some(held.iter().all(|hit| {
            rest.iter()
                .any(|other| other.marker == hit.marker && same(hit, other))
        }))
    }

    fn whole(&self, contig: usize, chosen: usize) -> Vec<usize> {
        let mut held = self.per_contig[contig]
            .iter()
            .filter(|hit| !hit.partial && self.set.sets.holds(chosen, hit.marker as usize))
            .filter(|hit| !self.is_multicopy(hit.marker))
            .map(|hit| hit.marker as usize)
            .collect::<Vec<_>>();
        held.dedup();
        held
    }

    fn counts(&self, contigs: &[usize]) -> Vec<Tally> {
        let mut counts = vec![Tally::default(); self.set.len()];
        for contig in contigs {
            for hit in &self.per_contig[*contig] {
                let tally = &mut counts[hit.marker as usize];
                tally.any += 1;
                tally.complete += u32::from(!hit.partial);
            }
        }
        if let Some(mask) = self.multicopy.get() {
            for (tally, multicopy) in counts.iter_mut().zip(mask) {
                if *multicopy {
                    tally.complete = tally.complete.min(1);
                }
            }
        }
        counts
    }

    fn is_multicopy(&self, marker: u16) -> bool {
        self.multicopy
            .get()
            .and_then(|mask| mask.get(marker as usize))
            .copied()
            .unwrap_or(false)
    }
}
