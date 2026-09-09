/// Which bar turned a re-clustering away. Counted per round so a refiner that splits nothing
/// says which test did it rather than leaving it to inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitRejection {
    SingleCluster,
    AllNoise,
    NotTighter,
    Shredded,
    Unimodal,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Rejections {
    pub too_few_contigs: usize,
    pub no_trigger: usize,
    pub no_clustering: usize,
    pub single_cluster: usize,
    pub all_noise: usize,
    pub not_tighter: usize,
    pub shredded: usize,
    pub unimodal: usize,
}

impl Rejections {
    pub fn merge(&mut self, other: &Self) {
        self.too_few_contigs += other.too_few_contigs;
        self.no_trigger += other.no_trigger;
        self.no_clustering += other.no_clustering;
        self.single_cluster += other.single_cluster;
        self.all_noise += other.all_noise;
        self.not_tighter += other.not_tighter;
        self.unimodal += other.unimodal;
        self.shredded += other.shredded;
    }

    pub fn record(&mut self, rejection: SplitRejection) {
        match rejection {
            SplitRejection::SingleCluster => self.single_cluster += 1,
            SplitRejection::AllNoise => self.all_noise += 1,
            SplitRejection::NotTighter => self.not_tighter += 1,
            SplitRejection::Shredded => self.shredded += 1,
            SplitRejection::Unimodal => self.unimodal += 1,
        }
    }
}

impl std::fmt::Display for Rejections {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "too few contigs {}, no trigger {}, no clustering {}, single cluster {}, \
             all noise {}, pieces not tighter {}, shredded {}, one mode {}",
            self.too_few_contigs,
            self.no_trigger,
            self.no_clustering,
            self.single_cluster,
            self.all_noise,
            self.not_tighter,
            self.shredded,
            self.unimodal
        )
    }
}

/// Why a bin was re-clustered, and for a tripped bin which test fired. Counted so a run that
/// splits everything says why, rather than leaving it to be inferred from the bins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trigger {
    Forced,
    Tripped { columns: [bool; 4], misplaced: bool },
    Peeled,
    Bisected,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct TriggerCounts {
    pub forced: usize,
    pub peeled: usize,
    pub bisected: usize,
    pub tripped: usize,
    pub columns: [usize; 4],
    pub misplaced: usize,
}

impl TriggerCounts {
    pub fn merge(&mut self, other: &Self) {
        self.forced += other.forced;
        self.peeled += other.peeled;
        self.bisected += other.bisected;
        self.tripped += other.tripped;
        for (total, add) in self.columns.iter_mut().zip(other.columns) {
            *total += add;
        }
        self.misplaced += other.misplaced;
    }

    pub fn record(&mut self, trigger: Trigger) {
        match trigger {
            Trigger::Forced => self.forced += 1,
            Trigger::Peeled => self.peeled += 1,
            Trigger::Bisected => self.bisected += 1,
            Trigger::Tripped {
                columns,
                misplaced: over_length,
            } => {
                self.tripped += 1;
                for (total, fired) in self.columns.iter_mut().zip(columns) {
                    *total += usize::from(fired);
                }
                self.misplaced += usize::from(over_length);
            }
        }
    }
}

impl std::fmt::Display for TriggerCounts {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "forced {}, peeled {}, bisected {}, tripped {} by metabat {}, \
             rho {}, euclidean {}, aggregate {}, misplaced length {}",
            self.forced,
            self.peeled,
            self.bisected,
            self.tripped,
            self.columns[0],
            self.columns[1],
            self.columns[2],
            self.columns[3],
            self.misplaced
        )
    }
}
