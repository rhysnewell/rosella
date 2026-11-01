use rayon::prelude::*;

use crate::embedding::knn::KnnGraph;
use crate::recover::fold::DepthFit;
use crate::refine::owners::{credit, owners};
use crate::refine::report_context::neighbour_weight;

const MOST: f64 = 0.5;

// A eukaryote barely carries the markers the arbiter weighs, so the partition cuts it freely.
// One round only, since a merged piece's few outside neighbours can pull in another genome.
pub fn rejoins(
    groups: &[&[usize]],
    eukaryotic: &[bool],
    depth: &DepthFit,
    knn: &KnnGraph,
    lengths: &[usize],
    shared: impl Fn(&[usize], &[usize]) -> Option<f64> + Sync,
) -> Vec<(usize, usize)> {
    let owner = owners(groups.iter().copied().enumerate());
    let pairs = (0..groups.len())
        .into_par_iter()
        .filter(|piece| eukaryotic[*piece])
        .filter_map(|piece| {
            let mut outside = Vec::new();
            for contig in groups[piece] {
                for (label, weight) in neighbour_weight(*contig, &owner, knn, lengths) {
                    if label != piece {
                        credit(&mut outside, label, weight);
                    }
                }
            }
            let total = outside.iter().map(|(_, weight)| weight).sum::<f64>();
            let (top, weight) = outside
                .into_iter()
                .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))?;
            (eukaryotic[top]
                && weight >= MOST * total
                && depth.fits(piece, top)
                && shared(groups[piece], groups[top]).is_none_or(|share| share <= MOST))
            .then_some((piece, top))
        })
        .collect::<Vec<_>>();
    roots(&pairs, groups, lengths)
}

// Each joined set goes to its largest piece, so no piece is moved into one that also moves.
fn roots(pairs: &[(usize, usize)], groups: &[&[usize]], lengths: &[usize]) -> Vec<(usize, usize)> {
    let mut parent = (0..groups.len()).collect::<Vec<_>>();
    fn find(parent: &mut [usize], mut at: usize) -> usize {
        while parent[at] != at {
            parent[at] = parent[parent[at]];
            at = parent[at];
        }
        at
    }
    let bases = |group: usize| {
        groups[group]
            .iter()
            .map(|contig| lengths[*contig])
            .sum::<usize>()
    };
    for (piece, top) in pairs {
        let (one, other) = (find(&mut parent, *piece), find(&mut parent, *top));
        if one == other {
            continue;
        }
        let keep = match bases(one).cmp(&bases(other)).then(other.cmp(&one)) {
            std::cmp::Ordering::Less => other,
            _ => one,
        };
        parent[one] = keep;
        parent[other] = keep;
    }
    let mut moved = (0..groups.len())
        .map(|group| (group, find(&mut parent, group)))
        .filter(|(group, root)| group != root)
        .collect::<Vec<_>>();
    moved.sort_unstable();
    moved
}
