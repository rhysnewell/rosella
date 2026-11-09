use std::collections::HashMap;

const LEAST: usize = 3;

// A pooled assembly holds one element once per sample. Its copies share their reads, so the
// mapper's split of those reads sets each copy's depth, and only ties and length name them.
pub struct Copies<'a> {
    pub lengths: &'a [usize],
    pub reads: &'a HashMap<usize, f64>,
    pub ties: &'a HashMap<(usize, usize), u64>,
    pub read_length: f64,
}

impl Copies<'_> {
    fn tied(&self, first: usize, second: usize, own: f64) -> bool {
        let shared = self
            .ties
            .get(&(first.min(second), first.max(second)))
            .copied()
            .unwrap_or(0);
        2.0 * shared as f64 >= own
    }

    fn copy(&self, first: usize, second: usize) -> bool {
        self.lengths[first].abs_diff(self.lengths[second]) as f64 <= self.read_length
            && self.tied(first, second, self.reads[&first].min(self.reads[&second]))
    }

    pub fn of(&self, members: &[usize]) -> Vec<Vec<usize>> {
        let mut root = (0..members.len()).collect::<Vec<_>>();
        fn find(root: &mut [usize], mut at: usize) -> usize {
            while root[at] != at {
                root[at] = root[root[at]];
                at = root[at];
            }
            at
        }
        for first in 0..members.len() {
            for second in first + 1..members.len() {
                if self.copy(members[first], members[second]) {
                    let (first, second) = (find(&mut root, first), find(&mut root, second));
                    root[first] = second;
                }
            }
        }
        let mut groups = HashMap::<usize, Vec<usize>>::new();
        for (at, contig) in members.iter().enumerate() {
            groups.entry(find(&mut root, at)).or_default().push(*contig);
        }
        let mut cores = groups
            .into_values()
            .filter(|group| {
                group.len() >= LEAST
                    && group.iter().enumerate().all(|(at, first)| {
                        group[at + 1..]
                            .iter()
                            .all(|second| self.copy(*first, *second))
                    })
            })
            .collect::<Vec<_>>();
        cores.sort_unstable();
        let sizes = cores.iter().map(Vec::len).collect::<Vec<_>>();
        for contig in members {
            if cores.iter().any(|core| core.contains(contig)) {
                continue;
            }
            let piece = cores.iter_mut().zip(&sizes).find(|(core, size)| {
                core[..**size].iter().all(|copy| {
                    self.lengths[*contig] < self.lengths[*copy]
                        && self.tied(*contig, *copy, self.reads[contig])
                })
            });
            if let Some((core, _)) = piece {
                core.push(*contig);
            }
        }
        cores
    }
}
