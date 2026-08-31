/// Which bar turned a re-clustering away. Counted per round so a refiner that splits nothing
/// says which test did it rather than leaving it to inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitRejection {
    SingleCluster,
    BelowTarget,
    AllNoise,
    NotTighter,
}

/// flight accepted a split on the density validity alone. `Strict` adds the noise cap and the
/// tightness test the port introduced, both of which can overrule that validity.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum SplitGate {
    #[default]
    Strict,
    Validity,
}

pub const SPLIT_GATE_NAMES: [&str; 2] = ["strict", "validity"];

impl SplitGate {
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "strict" => Some(Self::Strict),
            "validity" => Some(Self::Validity),
            _ => None,
        }
    }

    pub fn is_strict(&self) -> bool {
        *self == Self::Strict
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Rejections {
    pub too_few_contigs: usize,
    pub no_trigger: usize,
    pub no_clustering: usize,
    pub single_cluster: usize,
    pub below_target: usize,
    pub all_noise: usize,
    pub not_tighter: usize,
}

impl Rejections {
    pub fn merge(&mut self, other: &Self) {
        self.too_few_contigs += other.too_few_contigs;
        self.no_trigger += other.no_trigger;
        self.no_clustering += other.no_clustering;
        self.single_cluster += other.single_cluster;
        self.below_target += other.below_target;
        self.all_noise += other.all_noise;
        self.not_tighter += other.not_tighter;
    }

    pub fn record(&mut self, rejection: SplitRejection) {
        match rejection {
            SplitRejection::SingleCluster => self.single_cluster += 1,
            SplitRejection::BelowTarget => self.below_target += 1,
            SplitRejection::AllNoise => self.all_noise += 1,
            SplitRejection::NotTighter => self.not_tighter += 1,
        }
    }
}

impl std::fmt::Display for Rejections {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "too few contigs {}, no trigger {}, no clustering {}, single cluster {}, \
             below target {}, all noise {}, pieces not tighter {}",
            self.too_few_contigs,
            self.no_trigger,
            self.no_clustering,
            self.single_cluster,
            self.below_target,
            self.all_noise,
            self.not_tighter
        )
    }
}

/// Why a bin was re-clustered, and for a tripped bin which test fired. Counted so a run that
/// splits everything says why, rather than leaving it to be inferred from the bins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trigger {
    Forced,
    Tripped { columns: [bool; 4], misplaced: bool },
}

#[derive(Debug, Default, Clone, Copy)]
pub struct TriggerCounts {
    pub forced: usize,
    pub tripped: usize,
    pub columns: [usize; 4],
    pub misplaced: usize,
}

impl TriggerCounts {
    pub fn merge(&mut self, other: &Self) {
        self.forced += other.forced;
        self.tripped += other.tripped;
        for (total, add) in self.columns.iter_mut().zip(other.columns) {
            *total += add;
        }
        self.misplaced += other.misplaced;
    }

    pub fn record(&mut self, trigger: Trigger) {
        match trigger {
            Trigger::Forced => self.forced += 1,
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
            "forced {}, tripped {} by metabat {}, rho {}, euclidean {}, \
             aggregate {}, misplaced length {}",
            self.forced,
            self.tripped,
            self.columns[0],
            self.columns[1],
            self.columns[2],
            self.columns[3],
            self.misplaced
        )
    }
}
