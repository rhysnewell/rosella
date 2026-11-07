use std::sync::{
    Mutex,
    atomic::{AtomicU32, AtomicUsize, Ordering},
};

use rayon::prelude::*;

use super::NeighbourList;

// Forward and reverse neighbour lists, split by whether the edge is new since the last
// pass. Every bucket is capped, so one flat allocation of that stride holds the whole pass.
pub(super) struct Candidates {
    stride: usize,
    new: Vec<u32>,
    new_rank: Vec<u32>,
    new_len: Vec<u32>,
    old: Vec<u32>,
    old_rank: Vec<u32>,
    old_len: Vec<u32>,
    incoming: Incoming,
}

struct Bucket<'a> {
    values: &'a mut [u32],
    ranks: &'a mut [u32],
    length: &'a mut u32,
}

impl Bucket<'_> {
    // An equal rank never displaces the root, so each bucket hears its offers in row order.
    fn offer(&mut self, value: u32, rank: u32) {
        let (values, ranks) = (&mut *self.values, &mut *self.ranks);
        let stride = values.len();
        let length = *self.length as usize;
        if length < stride {
            *self.length += 1;
            values[length] = value;
            ranks[length] = rank;
            let mut at = length;
            while at > 0 && ranks[(at - 1) / 2] < ranks[at] {
                values.swap(at, (at - 1) / 2);
                ranks.swap(at, (at - 1) / 2);
                at = (at - 1) / 2;
            }
        } else if rank < ranks[0] {
            values[0] = value;
            ranks[0] = rank;
            let mut at = 0;
            loop {
                let largest = [2 * at + 1, 2 * at + 2]
                    .into_iter()
                    .filter(|child| *child < stride)
                    .fold(at, |best, child| match ranks[child] > ranks[best] {
                        true => child,
                        false => best,
                    });
                if largest == at {
                    break;
                }
                values.swap(at, largest);
                ranks.swap(at, largest);
                at = largest;
            }
        }
    }
}

// An edge holds its source above its new flag, so a sorted run reads in source order.
struct Incoming {
    starts: Vec<usize>,
    cursor: Vec<AtomicUsize>,
    edges: Vec<AtomicU32>,
}

impl Incoming {
    fn fill(&mut self, rows: &[Mutex<NeighbourList>]) {
        let n = rows.len();
        self.cursor
            .par_iter()
            .for_each(|count| count.store(0, Ordering::Relaxed));
        rows.par_iter().for_each(|row| {
            for to in &row.lock().unwrap().indices {
                if *to != u32::MAX {
                    self.cursor[*to as usize].fetch_add(1, Ordering::Relaxed);
                }
            }
        });
        self.starts[0] = 0;
        for at in 0..n {
            let count = self.cursor[at].load(Ordering::Relaxed);
            self.starts[at + 1] = self.starts[at] + count;
            self.cursor[at].store(self.starts[at], Ordering::Relaxed);
        }
        if self.edges.len() < self.starts[n] {
            self.edges.resize_with(self.starts[n], || AtomicU32::new(0));
        }
        rows.par_iter().enumerate().for_each(|(from, row)| {
            let list = row.lock().unwrap();
            for (to, fresh) in list.indices.iter().zip(&list.is_new) {
                if *to != u32::MAX {
                    let at = self.cursor[*to as usize].fetch_add(1, Ordering::Relaxed);
                    let edge = ((from as u32) << 1) | u32::from(*fresh);
                    self.edges[at].store(edge, Ordering::Relaxed);
                }
            }
        });
    }

    fn sorted_into(&self, row: usize, out: &mut Vec<u32>) {
        out.clear();
        out.extend(
            self.edges[self.starts[row]..self.starts[row + 1]]
                .iter()
                .map(|edge| edge.load(Ordering::Relaxed)),
        );
        out.sort_unstable();
    }
}

impl Candidates {
    pub(super) fn new(n: usize, stride: usize) -> Self {
        Self {
            stride,
            new: vec![u32::MAX; n * stride],
            new_rank: vec![u32::MAX; n * stride],
            new_len: vec![0; n],
            old: vec![u32::MAX; n * stride],
            old_rank: vec![u32::MAX; n * stride],
            old_len: vec![0; n],
            incoming: Incoming {
                starts: vec![0; n + 1],
                cursor: (0..n).map(|_| AtomicUsize::new(0)).collect(),
                edges: Vec::new(),
            },
        }
    }

    pub(super) fn new_of(&self, row: usize) -> &[u32] {
        &self.new[row * self.stride..row * self.stride + self.new_len[row] as usize]
    }

    pub(super) fn old_of(&self, row: usize) -> &[u32] {
        &self.old[row * self.stride..row * self.stride + self.old_len[row] as usize]
    }

    // Dong, Charikar & Li (2011) sample each bucket at random. Filled first come by row, late rows
    // starved and the search converged as far as the contig order let it.
    pub(super) fn build(&mut self, rows: &[Mutex<NeighbourList>], seed: u64) {
        self.incoming.fill(rows);
        let rank = |to: u32, from: usize| {
            (crate::digest::mix(seed ^ ((u64::from(to) << 32) | from as u64)) >> 32) as u32
        };
        let (stride, incoming) = (self.stride, &self.incoming);
        self.new
            .par_chunks_mut(stride)
            .zip(self.new_rank.par_chunks_mut(stride))
            .zip(self.new_len.par_iter_mut())
            .zip(
                self.old
                    .par_chunks_mut(stride)
                    .zip(self.old_rank.par_chunks_mut(stride))
                    .zip(self.old_len.par_iter_mut()),
            )
            .enumerate()
            .for_each_init(
                Vec::new,
                |sorted, (row, (((new, new_rank), new_len), ((old, old_rank), old_len)))| {
                    (*new_len, *old_len) = (0, 0);
                    let mut buckets = [
                        Bucket {
                            values: old,
                            ranks: old_rank,
                            length: old_len,
                        },
                        Bucket {
                            values: new,
                            ranks: new_rank,
                            length: new_len,
                        },
                    ];
                    incoming.sorted_into(row, sorted);
                    let split = sorted.partition_point(|edge| ((edge >> 1) as usize) < row);
                    let reverse = |buckets: &mut [Bucket; 2], edges: &[u32]| {
                        for edge in edges {
                            let from = (edge >> 1) as usize;
                            buckets[(edge & 1) as usize].offer(from as u32, rank(row as u32, from));
                        }
                    };
                    reverse(&mut buckets, &sorted[..split]);
                    {
                        let mut list = rows[row].lock().unwrap();
                        list.fresh.fill(false);
                        for (to, fresh) in list.indices.iter().zip(&list.is_new) {
                            if *to != u32::MAX {
                                buckets[usize::from(*fresh)].offer(*to, rank(*to, row));
                            }
                        }
                    }
                    reverse(&mut buckets, &sorted[split..]);
                },
            );
        self.dedup();
        rows.par_iter().enumerate().for_each(|(i, entry)| {
            let mut list = entry.lock().unwrap();
            let drawn = self.new_of(i);
            for slot in 0..list.indices.len() {
                if list.is_new[slot] && drawn.binary_search(&list.indices[slot]).is_ok() {
                    list.is_new[slot] = false;
                }
            }
        });
    }

    fn dedup(&mut self) {
        self.new
            .par_chunks_mut(self.stride)
            .zip(self.new_len.par_iter_mut())
            .zip(
                self.old
                    .par_chunks_mut(self.stride)
                    .zip(self.old_len.par_iter_mut()),
            )
            .for_each(|((new, new_len), (old, old_len))| {
                *new_len = unique(&mut new[..*new_len as usize], &[]);
                *old_len = unique(&mut old[..*old_len as usize], &new[..*new_len as usize]);
            });
    }
}

// A mutual edge reaches a bucket from both ends and a contig can sit in both buckets, so the join
// would measure one pair twice. The caps are already spent, so the lists come out the same.
fn unique(values: &mut [u32], elsewhere: &[u32]) -> u32 {
    values.sort_unstable();
    let mut kept = 0;
    for at in 0..values.len() {
        let value = values[at];
        if (kept == 0 || values[kept - 1] != value) && elsewhere.binary_search(&value).is_err() {
            values[kept] = value;
            kept += 1;
        }
    }
    kept as u32
}
