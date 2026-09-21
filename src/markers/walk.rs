use super::ContigMarkers;

/// A contig the shed evicted, with the carrier that already held its markers.
pub struct Shed {
    pub contig: usize,
    pub markers: usize,
    pub twin: Option<usize>,
    pub shared: usize,
}

impl ContigMarkers {
    /// A contig whose every whole marker copy the bin also holds on another contig cannot be
    /// carrying completeness, so it is either a second strain or a fragment of a neighbour.
    pub fn redundant(&self, contigs: &[usize], spacings: f64) -> Vec<usize> {
        self.walk(contigs, false, spacings)
            .into_iter()
            .map(|entry| entry.contig)
            .collect()
    }

    /// The twin search is off the shipped path, because naming the carrier that made a contig
    /// look redundant costs a pass over the bin for every eviction and only a probe reads it.
    pub fn redundant_traced(&self, contigs: &[usize], spacings: f64) -> Vec<Shed> {
        self.walk(contigs, true, spacings)
    }

    fn walk(&self, contigs: &[usize], trace: bool, spacings: f64) -> Vec<Shed> {
        let counts = self.counts(contigs);
        let Some(chosen) = self
            .set
            .sets
            .choose(&super::observed(&counts), self.bin_bp(contigs))
        else {
            return Vec::new();
        };
        let mut held = contigs.to_vec();
        let mut shed = Vec::new();
        while let Some(position) = self.passenger(&held, chosen, spacings) {
            let contig = held[position];
            let entry = match trace {
                true => self.trace(&held, contig, chosen),
                false => Shed {
                    contig,
                    markers: 0,
                    twin: None,
                    shared: 0,
                },
            };
            shed.push(entry);
            held.swap_remove(position);
        }
        shed.sort_unstable_by_key(|entry| entry.contig);
        shed
    }

    fn trace(&self, held: &[usize], contig: usize, chosen: usize) -> Shed {
        let mine = self.whole(contig, chosen);
        let mut best: Option<(usize, usize)> = None;
        for other in held.iter().filter(|other| **other != contig) {
            let shared = self
                .whole(*other, chosen)
                .iter()
                .filter(|marker| mine.contains(marker))
                .count();
            if shared > 0 && best.is_none_or(|(seen, _)| shared > seen) {
                best = Some((shared, *other));
            }
        }
        Shed {
            contig,
            markers: mine.len(),
            twin: best.map(|(_, other)| other),
            shared: best.map_or(0, |(shared, _)| shared),
        }
    }

    /// Carriers rather than copies, so the only contig holding a marker is never the one that
    /// leaves however many times it holds it. A contig long enough to be a stretch of the host
    /// genome is not a passenger, whatever its markers say, because the shed cannot tell a
    /// duplicated region from a second organism.
    fn passenger(&self, contigs: &[usize], chosen: usize, spacings: f64) -> Option<usize> {
        let bar = match spacings > 0.0 {
            true => spacings * self.set.sets.expected_bp(chosen)
                / self.set.sets.size(chosen).max(1) as f64,
            false => f64::INFINITY,
        };
        let mut carriers = vec![0u32; self.set.len()];
        for contig in contigs {
            for marker in self.whole(*contig, chosen) {
                carriers[marker] += 1;
            }
        }
        let mut best: Option<(usize, (usize, usize, usize))> = None;
        for (position, contig) in contigs.iter().enumerate() {
            if self.lengths.get(*contig).is_some_and(|bp| *bp as f64 >= bar) {
                continue;
            }
            let held = self.whole(*contig, chosen);
            if held.is_empty() || held.iter().any(|marker| carriers[*marker] < 2) {
                continue;
            }
            let key = (
                usize::MAX - held.len(),
                self.lengths.get(*contig).copied().unwrap_or_default(),
                *contig,
            );
            if best.is_none_or(|(_, seen)| key < seen) {
                best = Some((position, key));
            }
        }
        best.map(|(position, _)| position)
    }
}
