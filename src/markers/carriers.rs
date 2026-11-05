use super::ContigMarkers;

impl ContigMarkers {
    // Two contigs that each carry a whole copy of one single-copy marker hold two genomes.
    pub fn carriers(&self, contigs: &[usize]) -> Vec<Vec<u32>> {
        let (checkm_family, extra) = self.checkm_families();
        let mut families = vec![Vec::new(); self.set.len() + extra];
        for (at, contig) in contigs.iter().enumerate() {
            let hits = self.per_contig.get(*contig).map_or(&[][..], Vec::as_slice);
            let copies = self.checkm.get(*contig).map_or(&[][..], Vec::as_slice);
            let mut whole = hits
                .iter()
                .filter(|hit| !hit.partial && !self.set.eukaryotic(hit.marker))
                .map(|hit| hit.marker as usize)
                .chain(
                    copies
                        .iter()
                        .filter(|copy| !self.set.checkm.reads_gtdb(copy.set.into()))
                        .filter_map(|copy| checkm_family[usize::from(copy.model)])
                        .map(|family| self.set.len() + family),
                )
                .collect::<Vec<_>>();
            whole.sort_unstable();
            whole.dedup();
            for marker in whole {
                families[marker].push(at as u32);
            }
        }
        families.retain(|family| family.len() > 1);
        families
    }

    // A CheckM model searched under a GTDB name, or of a family GTDB holds, already reaches the
    // steps as that GTDB marker, so only the families GTDB lacks are added.
    fn checkm_families(&self) -> (Vec<Option<usize>>, usize) {
        let panel = &self.set.checkm;
        let mut index = std::collections::HashMap::<&str, usize>::new();
        let families = (0..panel.len() as u16)
            .map(|model| {
                let shared = panel
                    .searched_name(model)
                    .is_some_and(|name| self.set.id(name).is_some());
                let family = panel.family(model);
                (!shared && !self.set.holds_family(family)).then(|| {
                    let next = index.len();
                    *index.entry(family).or_insert(next)
                })
            })
            .collect();
        (families, index.len())
    }
}
