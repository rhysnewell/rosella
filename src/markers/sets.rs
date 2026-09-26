const STEPS: usize = 50;
const CLAMP: f64 = 1e-4;

#[derive(Default)]
pub struct Sets {
    names: Vec<String>,
    member_of: Vec<Vec<bool>>,
    sizes: Vec<usize>,
    present_log: Vec<Vec<f64>>,
    absent_log: Vec<Vec<f64>>,
    absent_total: Vec<Vec<f64>>,
    duplicate_weight: Vec<Vec<f64>>,
    max_bp: Vec<f64>,
    markers: usize,
}

impl Sets {
    pub fn new(
        names: Vec<String>,
        member_of: Vec<Vec<bool>>,
        rates: &[Vec<f64>],
        copies: Vec<Vec<f64>>,
    ) -> Self {
        let markers = member_of.first().map(Vec::len).unwrap_or_default();
        let sizes = member_of
            .iter()
            .map(|flags| flags.iter().filter(|held| **held).count())
            .collect();
        let mut present_log = Vec::with_capacity(names.len());
        let mut absent_log = Vec::with_capacity(names.len());
        let mut absent_total = Vec::with_capacity(names.len());
        for rate in rates {
            present_log.push(rate.iter().map(|u| clamp(*u).ln()).collect::<Vec<_>>());
            let mut rows = Vec::with_capacity(STEPS * markers);
            let mut totals = Vec::with_capacity(STEPS);
            for step in 1..=STEPS {
                let share = step as f64 / STEPS as f64;
                let start = rows.len();
                rows.extend(rate.iter().map(|u| (1.0 - clamp(share * u)).ln()));
                totals.push(rows[start..].iter().sum());
            }
            absent_log.push(rows);
            absent_total.push(totals);
        }
        Self {
            names,
            member_of,
            sizes,
            present_log,
            absent_log,
            absent_total,
            duplicate_weight: copies,
            max_bp: Vec::new(),
            markers,
        }
    }

    /// A second copy of a marker that is single copy in 86 per cent of genomes is weaker
    /// evidence of contamination than a second copy of one that is single copy in 99.
    pub fn duplicate_weight(&self, set: usize, marker: usize) -> f64 {
        self.duplicate_weight
            .get(set)
            .and_then(|held| held.get(marker))
            .copied()
            .unwrap_or(1.0)
    }

    pub fn with_bounds(mut self, max_bp: Vec<f64>) -> Self {
        self.max_bp = max_bp;
        self
    }

    fn fits(&self, set: usize, bin_bp: usize) -> bool {
        match self.max_bp.get(set).copied().unwrap_or_default() {
            bound if bound > 0.0 => bin_bp as f64 <= bound,
            _ => true,
        }
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn name(&self, set: usize) -> &str {
        self.names.get(set).map(String::as_str).unwrap_or_default()
    }

    pub fn holds(&self, set: usize, marker: usize) -> bool {
        self.member_of
            .get(set)
            .and_then(|flags| flags.get(marker))
            .copied()
            .unwrap_or_default()
    }

    pub fn size(&self, set: usize) -> usize {
        self.sizes.get(set).copied().unwrap_or_default()
    }

    pub fn widest(&self) -> Option<usize> {
        (0..self.len()).max_by_key(|set| self.sizes[*set])
    }

    /// A bin is a fraction of a genome, so each set is given its own completeness before the
    /// sets are compared. Without that the set with the lowest rates wins every sparse bin.
    /// A bin holds at least the genome it came from, so a set no genome that size belongs to
    /// is not a candidate however well its absences fit.
    pub fn choose(&self, present: &[u16], bin_bp: usize) -> Option<usize> {
        if self.markers == 0 || present.len() < crate::tuning::MARKERS_TO_CHOOSE_A_SET {
            return self.widest();
        }
        let mut best: Option<(usize, f64)> = None;
        for set in 0..self.len() {
            if self.sizes[set] == 0 || !self.fits(set, bin_bp) {
                continue;
            }
            let base = present
                .iter()
                .map(|marker| self.present_log[set][*marker as usize])
                .sum::<f64>();
            let mut top = f64::NEG_INFINITY;
            for step in 1..=STEPS {
                let row = &self.absent_log[set][(step - 1) * self.markers..step * self.markers];
                let taken = present
                    .iter()
                    .map(|marker| row[*marker as usize])
                    .sum::<f64>();
                let share = (step as f64 / STEPS as f64).ln();
                let total =
                    present.len() as f64 * share + base + self.absent_total[set][step - 1] - taken;
                if total > top {
                    top = total;
                }
            }
            if best.is_none_or(|(_, held)| top > held) {
                best = Some((set, top));
            }
        }
        best.map(|(set, _)| set).or_else(|| self.widest())
    }
}

fn clamp(rate: f64) -> f64 {
    rate.clamp(CLAMP, 1.0 - CLAMP)
}
