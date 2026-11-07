use std::cmp::Reverse;

use super::ContigMarkers;

impl ContigMarkers {
    // A contig whose every whole marker copy the bin also holds on another contig cannot be
    // carrying completeness, so it is either a second strain or a fragment of a neighbour.
    pub fn redundant(&self, contigs: &[usize]) -> Vec<usize> {
        let Some((chosen, counts)) = self.chosen(contigs) else {
            return Vec::new();
        };
        let mut held = contigs
            .iter()
            .map(|contig| (*contig, self.whole(*contig, chosen)))
            .collect::<Vec<_>>();
        if self.reads_checkm(contigs, self.read_on(&counts, chosen)) {
            for ((_, markers), models) in held.iter_mut().zip(self.checkm_carried(contigs, chosen))
            {
                markers.extend(models);
            }
        }
        let mut carriers = vec![0u32; self.set.len() + self.set.checkm.len()];
        for marker in held.iter().flat_map(|(_, markers)| markers) {
            carriers[*marker] += 1;
        }
        let mut shed = Vec::new();
        while let Some(position) = self.passenger(&held, &carriers) {
            let (contig, markers) = held.swap_remove(position);
            for marker in markers {
                carriers[marker] -= 1;
            }
            shed.push(contig);
        }
        shed.sort_unstable();
        shed
    }

    // Carriers rather than copies, so the only contig holding a marker is never the one that
    // leaves however many times it holds it.
    fn passenger(&self, held: &[(usize, Vec<usize>)], carriers: &[u32]) -> Option<usize> {
        held.iter()
            .enumerate()
            .filter(|(_, (_, markers))| {
                !markers.is_empty() && markers.iter().all(|marker| carriers[*marker] >= 2)
            })
            .min_by_key(|(_, (contig, markers))| {
                let length = self.lengths.get(*contig).copied().unwrap_or_default();
                (Reverse(markers.len()), length, *contig)
            })
            .map(|(position, _)| position)
    }
}
