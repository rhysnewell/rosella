use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::quality::Quality;

mod annotation;
mod annotator;
pub mod cache;
pub mod checkm;
mod community;
pub mod eukaryote;
mod fair;
pub mod fragments;
mod hit;
pub mod hmm_table;
mod multicopy;
pub mod panels;
mod readings;
pub mod replicon;
pub mod sets;
mod strain;
mod table;
mod walk;

pub use annotation::MarkerAnnotation;
pub use annotator::Annotator;
pub use fair::Fair;
pub use hit::{Hit, Place};
pub use readings::{CHECKM, COPY_BINS, GTDB, Reading};
pub use strain::{heterogeneity, identity, rebuild};
pub use table::MarkerSet;

pub(crate) const HMM_GZ: &[u8] = include_bytes!("../../data/gtdb_markers.hmm.gz");

const DOMAINS: [&str; 2] = ["bac", "ar"];

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
    checkm: Vec<Vec<checkm::Copies>>,
    shapes: Vec<replicon::Shape>,
    annotated: Vec<bool>,
    lengths: Vec<usize>,
    set: MarkerSet,
    multicopy: std::sync::OnceLock<Vec<bool>>,
    checkm_multicopy: std::sync::OnceLock<Vec<bool>>,
    community: std::sync::OnceLock<community::Community>,
    in_place: bool,
    fair: bool,
    sparse: Vec<bool>,
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
        self.score_by(contigs, self.fair)
    }

    fn points(&self, contigs: &[usize], weight: f64) -> Vec<(usize, f64)> {
        ContigMarkers::points(self, contigs, weight)
    }

    fn score_on(&self, contigs: &[usize], set: u16) -> Quality {
        self.score_on_by(contigs, set.into(), self.fair)
    }

    fn points_on(&self, contigs: &[usize], weight: f64, set: u16) -> Vec<(usize, f64)> {
        self.points_from(&self.counts(contigs), set.into(), weight)
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
            checkm_multicopy: std::sync::OnceLock::new(),
            community: std::sync::OnceLock::new(),
            in_place: false,
            fair: false,
            sparse: Vec::new(),
            per_contig,
            checkm: Vec::new(),
            shapes: Vec::new(),
            annotated: Vec::new(),
            lengths: Vec::new(),
            set,
        }
    }

    // Two contigs that each carry a whole copy of one single-copy marker hold two genomes.
    pub fn carriers(&self, contigs: &[usize]) -> Vec<Vec<u32>> {
        let mut families = vec![Vec::new(); self.set.len()];
        for (at, contig) in contigs.iter().enumerate() {
            let hits = self.per_contig.get(*contig).map_or(&[][..], Vec::as_slice);
            let mut whole = hits
                .iter()
                .filter(|hit| !hit.partial && !self.set.eukaryotic(hit.marker))
                .map(|hit| hit.marker as usize)
                .collect::<Vec<_>>();
            whole.dedup();
            for marker in whole {
                families[marker].push(at as u32);
            }
        }
        families.retain(|family| family.len() > 1);
        families
    }

    // Each genome carries one whole copy of a domain's single-copy markers, so the median
    // family's carrier count is how many genomes of that domain the assembly holds.
    pub fn genomes(&self) -> usize {
        let mut carriers = vec![0usize; self.set.len()];
        for hits in &self.per_contig {
            let mut whole = hits
                .iter()
                .filter(|hit| !hit.partial)
                .map(|hit| hit.marker as usize)
                .collect::<Vec<_>>();
            whole.dedup();
            for marker in whole {
                carriers[marker] += 1;
            }
        }
        let sets = &self.set.sets;
        (0..sets.len())
            .filter(|set| DOMAINS.contains(&sets.name(*set)))
            .map(|set| {
                let mut held = (0..carriers.len())
                    .filter(|marker| sets.holds(set, *marker))
                    .map(|marker| carriers[marker])
                    .collect::<Vec<_>>();
                held.sort_unstable();
                held.get(held.len() / 2).copied().unwrap_or_default()
            })
            .sum()
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

    pub fn with_checkm(mut self, checkm: Vec<Vec<checkm::Copies>>) -> Self {
        self.checkm = checkm;
        self
    }

    pub fn departures<'a>(
        &self,
        bins: impl IntoIterator<Item = &'a [usize]>,
        composition: impl Fn(&[usize]) -> Result<ndarray::Array2<f64>>,
        depths: ndarray::ArrayView2<f64>,
        depth_gated: bool,
    ) -> Result<replicon::Departures> {
        if self.shapes.len() != self.lengths.len() {
            return Ok(replicon::Departures::default());
        }
        let carries = (0..self.per_contig.len())
            .map(|contig| self.carries(contig))
            .collect::<Vec<_>>();
        replicon::departures(
            &self.shapes,
            &carries,
            &self.lengths,
            bins,
            composition,
            depths,
            depth_gated,
        )
    }

    fn chosen(&self, contigs: &[usize]) -> Option<(usize, Vec<Tally>)> {
        self.chosen_by(contigs, self.fair)
    }

    fn chosen_by(&self, contigs: &[usize], fair: bool) -> Option<(usize, Vec<Tally>)> {
        let counts = self.counted(contigs, fair);
        let chosen = self
            .set
            .sets
            .choose(&observed(&counts), self.bin_bp(contigs))?;
        Some((chosen, counts))
    }

    pub fn points(&self, contigs: &[usize], weight: f64) -> Vec<(usize, f64)> {
        self.points_by(contigs, weight, self.fair)
    }

    fn points_by(&self, contigs: &[usize], weight: f64, fair: bool) -> Vec<(usize, f64)> {
        let Some((chosen, counts)) = self.chosen_by(contigs, fair) else {
            return Vec::new();
        };
        self.points_from(&counts, chosen, weight)
    }

    fn points_from(&self, counts: &[Tally], set: usize, weight: f64) -> Vec<(usize, f64)> {
        counts
            .iter()
            .enumerate()
            .filter(|(marker, _)| self.set.sets.holds(set, *marker))
            .map(|(marker, tally)| {
                let extra = f64::from(tally.complete.saturating_sub(1));
                let present = f64::from(u8::from(tally.any >= 1));
                (marker, 100.0 * (present - weight * extra))
            })
            .collect()
    }

    // A eukaryotic panel hit says nothing of a prokaryote, so it never makes a carrier.
    pub fn carries(&self, contig: usize) -> bool {
        self.per_contig[contig]
            .iter()
            .any(|hit| !self.set.eukaryotic(hit.marker))
    }

    pub fn hit_count(&self, contigs: &[usize]) -> usize {
        contigs
            .iter()
            .flat_map(|contig| &self.per_contig[*contig])
            .filter(|hit| !self.set.eukaryotic(hit.marker))
            .count()
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
        self.counted(contigs, self.fair)
    }

    fn counted(&self, contigs: &[usize], fair: bool) -> Vec<Tally> {
        let mut counts = vec![Tally::default(); self.set.len()];
        for contig in contigs {
            let presence_only = self.presence_only(*contig);
            for hit in &self.per_contig[*contig] {
                let eukaryotic = self.set.eukaryotic(hit.marker);
                let tally = &mut counts[hit.marker as usize];
                tally.any += 1;
                tally.complete += u32::from(!hit.partial && (eukaryotic || !presence_only));
            }
        }
        if self.in_place {
            self.regroup(contigs, &mut counts, |hit| !self.set.eukaryotic(hit.marker));
        } else if fair {
            self.pair_cuts(contigs, &mut counts);
        }
        // Each exon of a gene split by introns can clear the cutoff alone, but only a second copy
        // covers the same stretch of the model.
        if self.set.has_panel() {
            self.regroup(contigs, &mut counts, |hit| {
                self.set.eukaryotic(hit.marker) && !hit.partial
            });
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

    // Bins from several binners differ by whole genomes cut into short contigs, whose doubled
    // markers are mostly cut genes, so a cut gene counts as a copy where it covers another's stretch.
    pub fn count_in_place(&mut self) {
        self.in_place = true;
    }

    fn regroup(&self, contigs: &[usize], counts: &mut [Tally], joins: impl Fn(&Hit) -> bool) {
        let mut held = std::collections::BTreeMap::<u16, Vec<&Hit>>::new();
        for hit in contigs.iter().flat_map(|contig| &self.per_contig[*contig]) {
            if counts[hit.marker as usize].any >= 2 && joins(hit) {
                held.entry(hit.marker).or_default().push(hit);
            }
        }
        for (marker, mut hits) in held {
            hits.sort_by_key(|hit| (hit.partial, hit.place.reach.model_from));
            let mut copies: Vec<Vec<&Hit>> = Vec::new();
            for hit in hits {
                match copies
                    .iter_mut()
                    .find(|copy| copy.iter().all(|other| !hit.same_part(other)))
                {
                    Some(copy) => copy.push(hit),
                    None => copies.push(vec![hit]),
                }
            }
            counts[marker as usize].complete = copies.len() as u32;
        }
    }

    fn is_multicopy(&self, marker: u16) -> bool {
        self.multicopy
            .get()
            .and_then(|mask| mask.get(marker as usize))
            .copied()
            .unwrap_or(false)
    }
}
