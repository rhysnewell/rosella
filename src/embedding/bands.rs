use ndarray::Array2;

/// Coverage cannot rank a pair it already places among `k` equals, so a gap inside that radius
/// is worth as little as a sample where neither contig is present.
#[derive(Debug)]
pub struct DepthBands {
    sorted: Vec<Vec<f64>>,
    k: usize,
    keeps_weight: bool,
}

impl DepthBands {
    pub fn new(coverage: &Array2<f64>, k: usize, keeps_weight: bool) -> Self {
        let n_samples = coverage.ncols() / 2;
        let sorted = (0..n_samples)
            .map(|sample| {
                let mut column = coverage
                    .column(sample * 2)
                    .iter()
                    .copied()
                    .collect::<Vec<_>>();
                column.sort_by(f64::total_cmp);
                column
            })
            .collect();
        Self {
            sorted,
            k: k.max(1),
            keeps_weight,
        }
    }

    pub fn keeps_weight(&self) -> bool {
        self.keeps_weight
    }

    pub fn radius(&self, sample: usize, depth: f64) -> f64 {
        let column = &self.sorted[sample];
        if column.len() <= self.k {
            return f64::INFINITY;
        }
        let at = column.partition_point(|value| *value < depth);
        let below = column[at.saturating_sub(self.k)];
        let above = column[(at + self.k).min(column.len() - 1)];
        (depth - below).min(above - depth)
    }
}
