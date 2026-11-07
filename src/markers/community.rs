use log::debug;

use super::eukaryote::LEAST_READINGS;
use super::multicopy::BINNED_COMPLETENESS;
use super::{ContigMarkers, Tally, checkm::Copies};
use crate::quality::{Bars, Quality, Scorer};

// Read once from the first ladder, so every stage after it judges a bin the same way.
#[derive(Clone, Copy, Debug)]
pub(super) struct Community {
    fence: Option<f64>,
    split: f32,
}

fn otsu(values: &mut [f64]) -> Option<f64> {
    if values.len() < LEAST_READINGS {
        return None;
    }
    values.sort_unstable_by(f64::total_cmp);
    let total = values.iter().sum::<f64>();
    let count = values.len() as f64;
    let (mut below, mut best) = (0.0, None::<(f64, usize)>);
    for at in 1..values.len() {
        below += values[at - 1];
        let (low, high) = (at as f64, count - at as f64);
        let gap = below / low - (total - below) / high;
        let between = low * high * gap * gap;
        if best.is_none_or(|(seen, _)| between > seen) {
            best = Some((between, at));
        }
    }
    best.map(|(_, at)| (values[at - 1] + values[at]) / 2.0)
}

impl ContigMarkers {
    pub fn learns_checkm(&self) -> bool {
        self.community.get().is_none()
    }

    // A paralog of a marker hits its model far under the genome's own copy, while a second
    // genome's copy scores close to it, so the extra copies of one community split in two.
    pub fn learn_checkm<'a>(&self, bins: impl IntoIterator<Item = &'a [usize]>, bars: Bars) {
        if !self.learns_checkm() {
            return;
        }
        let bins = bins.into_iter().collect::<Vec<_>>();
        let fence = self.eukaryote_fence(bins.iter().copied(), bars);
        let mut ratios = Vec::new();
        for contigs in &bins {
            let quality = self.score(contigs);
            if quality.completeness < BINNED_COMPLETENESS
                || self.reads_eukaryotic(contigs, quality, fence)
            {
                continue;
            }
            if let Some(at) = self.searched_lineage(quality.set.into()) {
                ratios.extend(self.extra_ratios(contigs, at).into_iter().map(f64::ln));
            }
        }
        let extras = ratios.len();
        let split = otsu(&mut ratios).map_or(0.0, |cut| cut.exp() as f32);
        debug!("CheckM paralog split: {split:.3} of the best copy ({extras} extra copies)");
        let _ = self.community.set(Community { fence, split });
    }

    pub(super) fn searched_lineage(&self, set: usize) -> Option<usize> {
        let panel = &self.set.checkm;
        panel
            .lineage(self.set.sets.name(set))
            .filter(|at| !panel.reads_gtdb(*at))
    }

    // GTDB under-reads contamination, so a prokaryote bin is judged on the higher of the two.
    pub(super) fn stricter(
        &self,
        contigs: &[usize],
        counts: &[Tally],
        quality: Quality,
    ) -> Quality {
        if !self.reads_checkm(contigs, quality) {
            return quality;
        }
        match self.checkm_read(contigs, counts, quality.set.into(), true) {
            Some(reading) if reading.contamination > quality.contamination => Quality {
                contamination: reading.contamination,
                ..quality
            },
            _ => quality,
        }
    }

    pub(super) fn reads_checkm(&self, contigs: &[usize], quality: Quality) -> bool {
        self.community
            .get()
            .is_some_and(|community| !self.reads_eukaryotic(contigs, quality, community.fence))
    }

    fn split(&self) -> f32 {
        self.community
            .get()
            .map_or(0.0, |community| community.split)
    }

    // A contig is one genome, so its copies of a model count once, at the best of them.
    pub(super) fn own_copies(&self, contig: usize, at: usize) -> impl Iterator<Item = &Copies> {
        let held = self.checkm.get(contig).map_or(&[][..], Vec::as_slice);
        held.iter().enumerate().filter_map(move |(position, copy)| {
            let first = position == 0
                || (held[position - 1].set, held[position - 1].model) != (copy.set, copy.model);
            (usize::from(copy.set) == at && first).then_some(copy)
        })
    }

    fn best_copies(&self, contigs: &[usize], at: usize) -> Vec<f32> {
        let mut best = vec![0.0f32; self.set.checkm.len()];
        for contig in contigs {
            for copy in self.own_copies(*contig, at) {
                let held = &mut best[usize::from(copy.model)];
                *held = held.max(copy.score);
            }
        }
        best
    }

    pub(super) fn each_counted(
        &self,
        contigs: &[usize],
        at: usize,
        mut visit: impl FnMut(usize, u16),
    ) {
        let best = self.best_copies(contigs, at);
        let split = self.split();
        for (position, contig) in contigs.iter().enumerate() {
            for copy in self.own_copies(*contig, at) {
                if copy.score >= split * best[usize::from(copy.model)] {
                    visit(position, copy.model);
                }
            }
        }
    }

    fn extra_ratios(&self, contigs: &[usize], at: usize) -> Vec<f64> {
        let best = self.best_copies(contigs, at);
        let mut seen = vec![false; best.len()];
        let mut ratios = Vec::new();
        for contig in contigs {
            for copy in self.own_copies(*contig, at) {
                let model = usize::from(copy.model);
                if copy.score >= best[model] && !seen[model] {
                    seen[model] = true;
                } else if best[model] > 0.0 {
                    ratios.push(f64::from(copy.score / best[model]));
                }
            }
        }
        ratios
    }
}
