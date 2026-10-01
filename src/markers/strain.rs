use std::collections::{BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader};

use anyhow::Result;
use log::info;
use rayon::prelude::*;

use super::{Annotator, ContigMarkers, Place, annotator, checkm};
use crate::external::hmmer_engine::HmmerEngine;
use crate::quality::orfs;

// CheckM1's --aai_strain default, compared strictly as CheckM1 does (Parks et al. 2015).
const STRAIN_IDENTITY: f64 = 0.9;

// Mirrors CheckM1's AAI, down to a backward scan that never reaches the first column.
pub fn identity(first: &[u8], second: &[u8]) -> f64 {
    let gap = |at: usize| first[at] == b'-' || second[at] == b'-';
    let length = first.len().min(second.len());
    let start = (0..length).find(|at| !gap(*at)).unwrap_or(length);
    let end = (1..length)
        .rev()
        .find(|at| !gap(*at))
        .map_or(length.min(1), |at| at + 1);
    let (mut mismatches, mut compared) = (0usize, 0usize);
    for at in start..end {
        if first[at] != second[at] {
            mismatches += 1;
            compared += 1;
        } else if first[at] != b'-' {
            compared += 1;
        }
    }
    match compared {
        0 => 0.0,
        _ => 1.0 - mismatches as f64 / compared as f64,
    }
}

pub fn heterogeneity(identities: &[f64]) -> f64 {
    match identities.len() {
        0 => 0.0,
        pairs => {
            let strains = identities
                .iter()
                .filter(|identity| **identity > STRAIN_IDENTITY)
                .count();
            100.0 * strains as f64 / pairs as f64
        }
    }
}

struct Wanted<'a> {
    slot: usize,
    at: usize,
    models: Vec<u16>,
    contigs: &'a [usize],
}

struct Held {
    slot: usize,
    profile: String,
    protein: String,
}

pub fn rebuild(sequence: &[u8], place: &Place) -> Option<String> {
    let begin = (place.gene_begin as usize).checked_sub(1)?;
    let coding = sequence.get(begin..place.gene_end as usize)?;
    Some(match place.reverse {
        true => orfs::translate_reverse(coding, !place.cut_right),
        false => orfs::translate(coding, !place.cut_left),
    })
}

impl Annotator {
    pub fn strain_heterogeneity(
        &self,
        markers: &ContigMarkers,
        bins: &[&[usize]],
        names: &[String],
    ) -> Result<Vec<Option<f64>>> {
        let mut out = vec![None; bins.len()];
        if !self.checkm {
            return Ok(out);
        }
        let (mut searched, mut read_off_gtdb) = (Vec::new(), Vec::new());
        for (slot, contigs) in bins.iter().enumerate() {
            let Some(doubled) = markers.doubled(contigs) else {
                continue;
            };
            out[slot] = Some(0.0);
            if doubled.models.is_empty() {
                continue;
            }
            let wanted = Wanted {
                slot,
                at: doubled.at,
                models: doubled.models,
                contigs,
            };
            match doubled.reads_gtdb {
                true => read_off_gtdb.push(wanted),
                false => searched.push(wanted),
            }
        }
        let mut copies = self.called_copies(markers, &searched, names)?;
        copies.extend(self.rebuilt_copies(markers, &read_off_gtdb, names)?);
        if copies.is_empty() {
            return Ok(out);
        }

        let _timer = crate::timing::scope("strain");
        let mut by_profile = HashMap::<String, Vec<(usize, String)>>::new();
        for held in copies {
            by_profile
                .entry(held.profile)
                .or_default()
                .push((held.slot, held.protein));
        }
        let profiles = profiles(&by_profile.keys().map(String::as_str).collect())?;
        let groups = by_profile.into_iter().collect::<Vec<_>>();
        let directory = tempfile::tempdir()?;
        let identities = groups
            .par_iter()
            .enumerate()
            .map(|(at, (name, copies))| match profiles.get(name) {
                Some(profile) => pairs(directory.path(), at, profile, copies),
                None => Ok(Vec::new()),
            })
            .collect::<Result<Vec<_>>>()?;
        let mut per_slot = HashMap::<usize, Vec<f64>>::new();
        for (slot, identity) in identities.into_iter().flatten() {
            per_slot.entry(slot).or_default().push(identity);
        }
        for (slot, identities) in per_slot {
            out[slot] = Some(heterogeneity(&identities));
        }
        Ok(out)
    }

    // A cached annotation holds no proteins, so the few contigs carrying a doubled marker are
    // called again.
    fn called_copies(
        &self,
        markers: &ContigMarkers,
        wanted: &[Wanted],
        names: &[String],
    ) -> Result<Vec<Held>> {
        let carries = |held: &Wanted, contig: usize| {
            markers.checkm[contig].iter().flatten().any(|entry| {
                usize::from(entry.set) == held.at && held.models.contains(&entry.model)
            })
        };
        let carriers = wanted
            .iter()
            .flat_map(|held| {
                held.contigs
                    .iter()
                    .copied()
                    .filter(move |contig| carries(held, *contig))
            })
            .collect::<BTreeSet<_>>();
        if carriers.is_empty() {
            return Ok(Vec::new());
        }
        let carriers = carriers.into_iter().collect::<Vec<_>>();
        info!(
            "Calling genes again on {} contigs that carry a doubled CheckM marker.",
            carriers.len()
        );
        let subset_names = carriers
            .iter()
            .map(|contig| names[*contig].clone())
            .collect::<Vec<_>>();
        let directory = tempfile::tempdir()?;
        let subset = directory.path().join("doubled.fna");
        annotator::write_contigs(&self.assembly, &subset_names, &subset)?;
        let annotator = Self {
            assembly: subset.to_string_lossy().into_owned(),
            cache: None,
            ..self.clone()
        };
        let mut called = annotator.annotate(0..usize::MAX)?.select(&subset_names)?;
        let mut homes = HashMap::<usize, Vec<&Wanted>>::new();
        for held in wanted {
            for contig in held.contigs {
                homes.entry(*contig).or_default().push(held);
            }
        }
        let panel = &markers.set.checkm;
        let mut found = Vec::new();
        for (copy, protein) in called.checkm_copies()? {
            for held in homes.get(&carriers[copy.contig]).into_iter().flatten() {
                let counted = usize::from(copy.set) == held.at && held.models.contains(&copy.model);
                if let Some(profile) = panel.searched_name(copy.model).filter(|_| counted) {
                    found.push(Held {
                        slot: held.slot,
                        profile: profile.to_string(),
                        protein: protein.clone(),
                    });
                }
            }
        }
        Ok(found)
    }

    // A GTDB hit keeps its gene's place, so its protein is translated again rather than called.
    fn rebuilt_copies(
        &self,
        markers: &ContigMarkers,
        wanted: &[Wanted],
        names: &[String],
    ) -> Result<Vec<Held>> {
        let mut places = Vec::new();
        for held in wanted {
            for contig in held.contigs {
                let hits = markers.per_contig[*contig]
                    .iter()
                    .filter(|hit| held.models.contains(&hit.marker));
                places.extend(hits.map(|hit| (held.slot, *contig, hit.marker, hit.place)));
            }
        }
        if places.is_empty() {
            return Ok(Vec::new());
        }
        let carriers = places
            .iter()
            .map(|(_, contig, _, _)| names[*contig].as_str())
            .collect::<HashSet<_>>();
        let mut sequences = HashMap::new();
        crate::kmers::each_named(
            &self.assembly,
            |name| carriers.contains(name),
            |name, record| {
                sequences.insert(name.to_string(), record.seq().into_owned());
                Ok(())
            },
        )?;
        Ok(places
            .into_iter()
            .filter_map(|(slot, contig, marker, place)| {
                Some(Held {
                    slot,
                    profile: markers.set.name(marker).to_string(),
                    protein: rebuild(sequences.get(&names[contig])?, &place)?,
                })
            })
            .collect())
    }
}

fn pairs(
    directory: &std::path::Path,
    at: usize,
    profile: &str,
    copies: &[(usize, String)],
) -> Result<Vec<(usize, f64)>> {
    let mut per_slot = HashMap::<usize, usize>::new();
    for (slot, _) in copies {
        *per_slot.entry(*slot).or_default() += 1;
    }
    if per_slot.values().all(|count| *count < 2) {
        return Ok(Vec::new());
    }
    let hmm = directory.join(format!("{at}.hmm"));
    let sequences = directory.join(format!("{at}.faa"));
    std::fs::write(&hmm, profile)?;
    let fasta = copies
        .iter()
        .enumerate()
        .map(|(at, (_, sequence))| format!(">{at}\n{sequence}\n"))
        .collect::<String>();
    std::fs::write(&sequences, fasta)?;
    let rows = masked(&HmmerEngine::align(&hmm, &sequences)?, copies.len());
    let mut found = Vec::new();
    for (first, (slot, _)) in copies.iter().enumerate() {
        for (second, (other, _)) in copies.iter().enumerate().skip(first + 1) {
            if slot == other {
                found.push((*slot, identity(&rows[first], &rows[second])));
            }
        }
    }
    Ok(found)
}

// Dropping A2M's insert columns leaves the match columns CheckM1's reference-line mask keeps.
fn masked(alignment: &str, count: usize) -> Vec<Vec<u8>> {
    let mut rows = vec![Vec::new(); count];
    let mut at = None;
    for line in alignment.lines() {
        if let Some(header) = line.strip_prefix('>') {
            at = header
                .split_whitespace()
                .next()
                .and_then(|id| id.parse::<usize>().ok());
            continue;
        }
        if let Some(row) = at.and_then(|at| rows.get_mut(at)) {
            row.extend(
                line.bytes()
                    .filter(|residue| residue.is_ascii_uppercase() || *residue == b'-'),
            );
        }
    }
    rows
}

fn profiles(wanted: &HashSet<&str>) -> Result<HashMap<String, String>> {
    let mut found = HashMap::new();
    for compressed in [super::HMM_GZ, checkm::HMM_GZ] {
        let reader = BufReader::new(flate2::read::GzDecoder::new(compressed));
        let mut block = String::new();
        let mut name = None;
        for line in reader.lines() {
            let line = line?;
            if let Some(rest) = line.strip_prefix("NAME") {
                name = Some(rest.trim().to_string());
            }
            block.push_str(&line);
            block.push('\n');
            if line == "//" {
                if let Some(held) = name.take().filter(|held| wanted.contains(held.as_str())) {
                    found.insert(held, std::mem::take(&mut block));
                }
                block.clear();
            }
        }
    }
    Ok(found)
}
