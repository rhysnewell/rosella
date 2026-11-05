use super::ContigMarkers;

impl ContigMarkers {
    // Two contigs that each carry a whole copy of one single-copy marker hold two genomes.
    pub fn carriers(&self, contigs: &[usize]) -> Vec<Vec<u32>> {
        let (checkm_family, extra) = self.checkm_families();
        let mut families = vec![Vec::new(); self.set.len() + extra];
        for (at, contig) in contigs.iter().enumerate() {
            let hits = self.per_contig.get(*contig).map_or(&[][..], Vec::as_slice);
            let copies = self.checkm.get(*contig).map_or(&[][..], Vec::as_slice);
            let lineage = self.carrier_lineage(*contig);
            let mut whole = hits
                .iter()
                .filter(|hit| !hit.partial && !self.set.eukaryotic(hit.marker))
                .map(|hit| hit.marker as usize)
                .chain(
                    copies
                        .iter()
                        .filter(|copy| !copy.partial && Some(usize::from(copy.set)) == lineage)
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

    // CheckM reads only the set of a genome's own lineage (Parks et al. 2015), and a domain that
    // is single copy in archaea can sit in several genes of one bacterium.
    fn carrier_lineage(&self, contig: usize) -> Option<usize> {
        let mut present = self
            .per_contig
            .get(contig)?
            .iter()
            .map(|hit| hit.marker)
            .collect::<Vec<_>>();
        present.dedup();
        let set = self.set.sets.place(&present, self.bin_bp(&[contig]))?;
        self.searched_lineage(set)
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
