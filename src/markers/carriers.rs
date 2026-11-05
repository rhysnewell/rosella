use super::ContigMarkers;

impl ContigMarkers {
    // A CheckM family GTDB lacks joins only while it sits on no more contigs than one GTDB holds,
    // since a domain one genome carries in several genes would part that genome from itself.
    pub fn carriers(&self, contigs: &[usize]) -> Vec<Vec<u32>> {
        let (family_of, held) = self.checkm_families();
        let mut families = vec![Vec::new(); self.set.len()];
        let mut checkm = vec![Vec::new(); held.len()];
        for (at, contig) in contigs.iter().enumerate() {
            let hits = self.per_contig.get(*contig).map_or(&[][..], Vec::as_slice);
            let whole = hits
                .iter()
                .filter(|hit| !hit.partial && !self.set.eukaryotic(hit.marker))
                .map(|hit| hit.marker as usize);
            let copies = self.checkm.get(*contig).map_or(&[][..], Vec::as_slice);
            let found = copies
                .iter()
                .filter(|copy| !self.set.checkm.reads_gtdb(copy.set.into()))
                .map(|copy| family_of[usize::from(copy.model)]);
            push_once(&mut families, whole, at);
            push_once(&mut checkm, found, at);
        }
        let ceiling = checkm
            .iter()
            .zip(&held)
            .filter_map(|(carriers, held)| held.then_some(carriers.len()))
            .max()
            .unwrap_or(0);
        families.extend(
            checkm
                .into_iter()
                .zip(held)
                .filter(|(carriers, held)| !held && carriers.len() <= ceiling)
                .map(|(carriers, _)| carriers),
        );
        families.retain(|family| family.len() > 1);
        families
    }

    fn checkm_families(&self) -> (Vec<usize>, Vec<bool>) {
        let panel = &self.set.checkm;
        let mut index = std::collections::HashMap::<&str, usize>::new();
        let mut held = Vec::new();
        let family_of = (0..panel.len() as u16)
            .map(|model| {
                let family = panel.family(model);
                let at = *index.entry(family).or_insert_with(|| {
                    held.push(false);
                    held.len() - 1
                });
                held[at] |= self.set.holds_family(family)
                    || panel
                        .searched_name(model)
                        .is_some_and(|name| self.set.id(name).is_some());
                at
            })
            .collect();
        (family_of, held)
    }
}

fn push_once(families: &mut [Vec<u32>], found: impl Iterator<Item = usize>, at: usize) {
    let mut found = found.collect::<Vec<_>>();
    found.sort_unstable();
    found.dedup();
    for family in found {
        families[family].push(at as u32);
    }
}
