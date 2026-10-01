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
    pub fn redundant(&self, contigs: &[usize]) -> Vec<usize> {
        self.walk(contigs, false)
            .into_iter()
            .map(|entry| entry.contig)
            .collect()
    }

    /// The twin search is off the shipped path, because naming the carrier that made a contig
    /// look redundant costs a pass over the bin for every eviction and only a probe reads it.
    pub fn redundant_traced(&self, contigs: &[usize]) -> Vec<Shed> {
        self.walk(contigs, true)
    }

    fn walk(&self, contigs: &[usize], trace: bool) -> Vec<Shed> {
        let Some((chosen, _)) = self.chosen(contigs) else {
            return Vec::new();
        };
        let mut held = contigs
            .iter()
            .map(|contig| (*contig, self.whole(*contig, chosen)))
            .collect::<Vec<_>>();
        let mut carriers = vec![0u32; self.set.len()];
        for marker in held.iter().flat_map(|(_, markers)| markers) {
            carriers[*marker] += 1;
        }
        let mut shed = Vec::new();
        while let Some(position) = self.passenger(&held, &carriers) {
            let (contig, markers) = &held[position];
            shed.push(match trace {
                true => twin_of(&held, *contig, markers),
                false => Shed {
                    contig: *contig,
                    markers: 0,
                    twin: None,
                    shared: 0,
                },
            });
            for marker in held.swap_remove(position).1 {
                carriers[marker] -= 1;
            }
        }
        shed.sort_unstable_by_key(|entry| entry.contig);
        shed
    }

    /// Carriers rather than copies, so the only contig holding a marker is never the one that
    /// leaves however many times it holds it.
    fn passenger(&self, held: &[(usize, Vec<usize>)], carriers: &[u32]) -> Option<usize> {
        let mut best: Option<(usize, (usize, usize, usize))> = None;
        for (position, (contig, markers)) in held.iter().enumerate() {
            if markers.is_empty() || markers.iter().any(|marker| carriers[*marker] < 2) {
                continue;
            }
            let key = (
                usize::MAX - markers.len(),
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

fn twin_of(held: &[(usize, Vec<usize>)], contig: usize, mine: &[usize]) -> Shed {
    let mut best: Option<(usize, usize)> = None;
    for (other, theirs) in held.iter().filter(|(other, _)| *other != contig) {
        let shared = theirs.iter().filter(|marker| mine.contains(marker)).count();
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
