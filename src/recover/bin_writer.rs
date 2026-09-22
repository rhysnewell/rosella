use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet, hash_map::Entry},
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    path,
};

use anyhow::{Context, Result};
use log::{debug, info, warn};
use needletail::{
    parse_fastx_file,
    parser::{LineEnding, write_fasta},
};
use rayon::slice::ParallelSliceMut;

use crate::recover::recover_engine::{RecoverEngine, UNBINNED};

pub const ELEMENT_PREFIX: &str = "viral_";

impl RecoverEngine {
    pub(crate) fn get_cluster_result(
        &self,
        cluster_map: HashMap<usize, HashSet<usize>>,
        outliers: HashSet<usize>,
    ) -> Vec<ClusterResult> {
        let elements = self.quality.small_elements();
        let mut cluster_results = Vec::with_capacity(self.n_contigs);
        for (cluster_label, contig_indices) in cluster_map.into_iter() {
            let bin_size = contig_indices
                .iter()
                .filter(|i| !elements.contains(*i))
                .map(|i| self.coverage_table.contig_lengths[*i])
                .sum::<usize>();
            let cluster_label = if bin_size < self.min_bin_size {
                None
            } else {
                Some(cluster_label)
            };
            for contig_index in contig_indices.iter() {
                cluster_results.push(if elements.contains(contig_index) {
                    ClusterResult::element(*contig_index)
                } else {
                    ClusterResult::new(*contig_index, cluster_label)
                });
            }
        }
        for outlier in outliers {
            cluster_results.push(if elements.contains(&outlier) {
                ClusterResult::element(outlier)
            } else {
                ClusterResult::new(outlier, None)
            });
        }
        cluster_results.par_sort_unstable();

        debug!(
            "Cluster results: {:?}",
            &cluster_results[..cluster_results.len().min(10)]
        );
        cluster_results
    }

    /// Take the cluster results and collect the contigs into bins.
    ///
    /// Keyed on contig name rather than position. The clustering indexes the coverage
    /// table, which the length filter has already shortened, so walking the assembly and
    /// counting sends every contig after the first short one to the wrong bin.
    pub(crate) fn write_clusters(&self, cluster_results: &[ClusterResult]) -> Result<()> {
        let labels = cluster_results
            .iter()
            .map(|result| {
                (
                    self.coverage_table.contig_names[result.contig_index].as_str(),
                    (result.cluster_label, result.element),
                )
            })
            .collect::<HashMap<_, _>>();

        let mut reader = parse_fastx_file(path::Path::new(&self.assembly))?;
        let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();
        let mut single_contig_bin_id = 0;
        let mut element_bin_id = 0;
        let mut unrecognised = 0;
        let mut read = 0;
        let mut written = 0;
        let progress = crate::progress::spinning(crate::progress::Stage::WritingBins);

        while let Some(record) = reader.next() {
            let seqrec = record?;
            read += 1;
            let contig_name = crate::contig_id(seqrec.id())?.to_string();
            let contig_length = seqrec.seq().len();

            let cluster_label = if contig_length < self.min_contig_size {
                self.leftover_label(contig_length, self.min_bin_size, &mut single_contig_bin_id)
            } else {
                match labels.get(contig_name.as_str()) {
                    Some((_, true)) => {
                        element_bin_id += 1;
                        format!("{ELEMENT_PREFIX}{element_bin_id}")
                    }
                    Some((Some(cluster_label), _)) => format!("{cluster_label}"),
                    Some((None, _)) => self.leftover_label(
                        contig_length,
                        self.min_bin_size,
                        &mut single_contig_bin_id,
                    ),
                    None => {
                        unrecognised += 1;
                        self.leftover_label(
                            contig_length,
                            self.min_bin_size,
                            &mut single_contig_bin_id,
                        )
                    }
                }
            };

            let writer = match writers.entry(cluster_label) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    let bin_path = path::Path::new(&self.output_directory).join(format!(
                        "rosella_bin_{}.{}",
                        entry.key(),
                        crate::defaults::FASTA_EXTENSION
                    ));
                    let file = OpenOptions::new()
                        .append(true)
                        .create(true)
                        .open(bin_path)?;
                    entry.insert(BufWriter::new(file))
                }
            };
            write_fasta(seqrec.id(), &seqrec.seq(), writer, LineEnding::Unix)?;
            written += 1;
            progress.set_message(format!("{written} contigs, {} bins", writers.len()));
        }
        progress.finish_and_clear();
        let n_bins = writers.len();

        // Dropping a BufWriter flushes it and throws the error away, so a full disk or a
        // broken pipe would truncate a bin silently.
        for (label, mut writer) in writers.drain() {
            writer
                .flush()
                .with_context(|| format!("flushing rosella_bin_{label}"))?;
        }

        if unrecognised > 0 {
            warn!(
                "{} contigs in the assembly were not in the coverage table and went unbinned",
                unrecognised
            );
        }
        if written != read {
            bail!(
                "{} of {} assembly contigs were written. Every contig belongs in a bin, in \
                 rosella_bin_unbinned or in rosella_bin_small_unbinned, so a shortfall means \
                 contigs were dropped",
                written,
                read
            );
        }
        info!(
            "Wrote {written} contigs into {n_bins} bins in {}.",
            self.output_directory
        );

        Ok(())
    }

    /// Where a contig goes when it has no cluster of its own: a bin by itself if it is
    /// long enough to be worth reporting, otherwise the unbinned pile.
    fn leftover_label(&self, contig_length: usize, floor: usize, next_id: &mut usize) -> String {
        if contig_length >= floor {
            *next_id += 1;
            return format!("single_contig_{next_id}");
        }
        if contig_length < self.min_contig_size {
            return format!("small_{UNBINNED}");
        }
        UNBINNED.to_string()
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct ClusterResult {
    pub(crate) contig_index: usize,
    pub(crate) cluster_label: Option<usize>,
    pub(crate) element: bool,
}

impl ClusterResult {
    pub fn new(contig_index: usize, cluster_label: Option<usize>) -> Self {
        Self {
            contig_index,
            cluster_label,
            element: false,
        }
    }

    pub fn element(contig_index: usize) -> Self {
        Self {
            contig_index,
            cluster_label: None,
            element: true,
        }
    }
}

/// Contig index first, then the label, so the unstable parallel sort at `write_clusters` has
/// no tie for the work stealing split to pick.
impl Ord for ClusterResult {
    fn cmp(&self, other: &Self) -> Ordering {
        self.contig_index
            .cmp(&other.contig_index)
            .then(self.cluster_label.cmp(&other.cluster_label))
    }
}

impl PartialOrd for ClusterResult {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
