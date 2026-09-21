/// The pool reads this on its first pass, before any stage takes a bin apart, so a later stage
/// can ask what kind of assembly it is on without a second look at the data.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Finished {
    held_back: usize,
    handed: usize,
}

impl Finished {
    pub fn new(held_back: usize, handed: usize) -> Self {
        Self { held_back, handed }
    }

    pub fn share(self) -> f64 {
        if self.handed == 0 {
            0.0
        } else {
            self.held_back as f64 / self.handed as f64
        }
    }

    /// Under a majority, sitting below the bars says nothing about a bin, since most of them do.
    /// Over one it marks a genome as incomplete rather than a pair as fused, and taking it apart
    /// then spends sequence looking for a second genome that was never in it.
    pub fn mostly(self) -> bool {
        self.handed > 0 && 2 * self.held_back >= self.handed
    }
}
