/// Which bar turned a re-clustering away. Counted per round so a refiner that splits nothing
/// says which test did it rather than leaving it to inference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitRejection {
    SingleCluster,
    BelowTarget,
    AllNoise,
    NoBinOverFloor,
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
    pub already_clean: usize,
    pub no_clustering: usize,
    pub single_cluster: usize,
    pub below_target: usize,
    pub all_noise: usize,
    pub no_bin_over_floor: usize,
    pub not_tighter: usize,
}

impl Rejections {
    pub fn record(&mut self, rejection: SplitRejection) {
        match rejection {
            SplitRejection::SingleCluster => self.single_cluster += 1,
            SplitRejection::BelowTarget => self.below_target += 1,
            SplitRejection::AllNoise => self.all_noise += 1,
            SplitRejection::NoBinOverFloor => self.no_bin_over_floor += 1,
            SplitRejection::NotTighter => self.not_tighter += 1,
        }
    }
}

impl std::fmt::Display for Rejections {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "too few contigs {}, already clean {}, no clustering {}, single cluster {}, \
             below target {}, all noise {}, no bin over the floor {}, pieces not tighter {}",
            self.too_few_contigs,
            self.already_clean,
            self.no_clustering,
            self.single_cluster,
            self.below_target,
            self.all_noise,
            self.no_bin_over_floor,
            self.not_tighter
        )
    }
}
