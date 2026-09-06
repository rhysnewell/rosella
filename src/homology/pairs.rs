#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pair {
    pub one: usize,
    pub other: usize,
    pub identity: f64,
    pub aligned_one: f64,
    pub aligned_other: f64,
}

impl Pair {
    /// Both ends have to be covered. A repeat inside one genome aligns over all of the short
    /// contig and almost none of the long one, so the smaller fraction is the honest number.
    pub fn aligned_fraction(&self) -> f64 {
        self.aligned_one.min(self.aligned_other)
    }
}
