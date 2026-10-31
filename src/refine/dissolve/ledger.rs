use crate::embedding::features::ContigFeatures;
use crate::refine::finished::Finished;

#[derive(Debug, Default, Clone, Copy)]
pub struct Mass {
    pub contigs: usize,
    pub bp: usize,
}

impl Mass {
    pub fn of<'a>(features: &ContigFeatures, contigs: impl IntoIterator<Item = &'a usize>) -> Self {
        contigs
            .into_iter()
            .fold(Self::default(), |mass, contig| Self {
                contigs: mass.contigs + 1,
                bp: mass.bp + features.length(*contig),
            })
    }

    pub fn add(&mut self, other: Self) {
        self.contigs += other.contigs;
        self.bp += other.bp;
    }
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Handed {
    pub held_back: usize,
    pub small: usize,
    pub clean: usize,
    pub bp: usize,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Search {
    pub rounds: usize,
    pub passes: usize,
    pub rung: usize,
    pub proposed: usize,
    pub composition: usize,
    pub linkage: usize,
    pub noise: usize,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Refused {
    pub small: usize,
    pub incomplete: usize,
    pub contaminated: usize,
    pub consumed: usize,
    pub worse: usize,
    pub carved: usize,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Claims {
    pub promoted: usize,
    pub deferred: usize,
    pub drained: usize,
    pub restored: usize,
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Placed {
    pub adopted: Mass,
    pub folded: Mass,
    pub returned: Mass,
    pub emptied: usize,
    pub left: Mass,
}

// What the pool took, what it refused and where the refusals went, in contigs and bases.
// The bin count alone cannot say whether an arm has headroom left.
#[derive(Debug, Default, Clone, Copy)]
pub struct DissolveLedger {
    pub handed: Handed,
    pub pool: Mass,
    pub search: Search,
    pub refused: Refused,
    pub claims: Claims,
    pub placed: Placed,
}

impl DissolveLedger {
    pub fn finished(&self) -> Finished {
        let handed = self.handed;
        Finished::new(
            handed.held_back,
            handed.held_back + handed.small + handed.clean,
        )
    }
}

impl std::fmt::Display for DissolveLedger {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            handed,
            pool,
            search,
            refused,
            claims,
            placed,
        } = self;
        write!(
            formatter,
            "held back {} bins already over the bars; dissolved {} small and {} \
             clean bins holding {} bp; pool {} contigs {} bp; {} rounds over {} passes ending at \
             rung {}; proposed {} clusters, refused {} small, {} incomplete, {} contaminated, {} \
             eaten by a rival, {} no better and {} carved out of a bin with no duplication to \
             explain them, {} noise; {} of the proposals came only from composition and the \
             merge order offered {}; promoted {} bins adopting {} contigs {} bp; folded {} \
             contigs {} bp back into the claims that took their bins; returned {} contigs {} bp, \
             emptied {} bins; left {} contigs {} bp unbinned; restored {} bins the pool broke \
             into nothing; held {} loose candidates back and drained {} of them",
            handed.held_back,
            handed.small,
            handed.clean,
            handed.bp,
            pool.contigs,
            pool.bp,
            search.rounds,
            search.passes,
            search.rung,
            search.proposed,
            refused.small,
            refused.incomplete,
            refused.contaminated,
            refused.consumed,
            refused.worse,
            refused.carved,
            search.noise,
            search.composition,
            search.linkage,
            claims.promoted,
            placed.adopted.contigs,
            placed.adopted.bp,
            placed.folded.contigs,
            placed.folded.bp,
            placed.returned.contigs,
            placed.returned.bp,
            placed.emptied,
            placed.left.contigs,
            placed.left.bp,
            claims.restored,
            claims.deferred,
            claims.drained
        )
    }
}
