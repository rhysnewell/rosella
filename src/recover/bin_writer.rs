use std::{
    cmp::Ordering,
    collections::{HashMap, HashSet, hash_map::Entry},
    fs::{File, OpenOptions},
    io::BufWriter,
    path,
};

use anyhow::Result;
use log::{debug, warn};
use needletail::{
    parse_fastx_file,
    parser::{LineEnding, write_fasta},
};
use rayon::slice::ParallelSliceMut;

use crate::recover::recover_engine::{
    RECOVER_FASTA_EXTENSION, REFINING_BIN_SIZE, RecoverEngine, UNBINNED,
};

impl RecoverEngine {
    pub(crate) fn get_cluster_result(
        &self,
        cluster_map: HashMap<usize, HashSet<usize>>,
        outliers: HashSet<usize>,
        contig_index_map: Option<HashMap<usize, usize>>, // a map containing the key as the row index in the embeddings array and the value as the original contig index
                                                         // used when a cluster has been subset and re-embedded
    ) -> Vec<ClusterResult> {
        let mut cluster_results = Vec::with_capacity(self.n_contigs);
        for (cluster_label, mut contig_indices) in cluster_map.into_iter() {
            contig_indices = match &contig_index_map {
                Some(index_map) => contig_indices
                    .into_iter()
                    .filter_map(|idx| index_map.get(&idx).map(|idx| *idx))
                    .collect::<HashSet<_>>(),
                None => contig_indices,
            };
            // check the size of the cluster and if it is too small, set the cluster label to None
            let bin_size = contig_indices
                .iter()
                .map(|i| self.coverage_table.contig_lengths[*i])
                .sum::<usize>();
            let cluster_label = if bin_size < self.min_bin_size {
                None
            } else {
                Some(cluster_label)
            };
            for contig_index in contig_indices.iter() {
                cluster_results.push(ClusterResult::new(*contig_index, cluster_label));
            }
        }
        for outlier in outliers {
            cluster_results.push(ClusterResult::new(outlier, None));
        }
        cluster_results.par_sort_unstable();

        debug!("Cluster results: {:?}", &cluster_results[0..10]);
        cluster_results
    }

    /// Take the cluster results and collect the contigs into bins.
    ///
    /// Keyed on contig name rather than position. The clustering indexes the coverage
    /// table, which the length filter has already shortened, so walking the assembly and
    /// counting sends every contig after the first short one to the wrong bin.
    pub(crate) fn write_clusters(
        &self,
        cluster_results: Vec<ClusterResult>,
        to_refine: bool,
    ) -> Result<()> {
        let labels = cluster_results
            .iter()
            .map(|result| {
                (
                    self.coverage_table.contig_names[result.contig_index].as_str(),
                    result.cluster_label,
                )
            })
            .collect::<HashMap<_, _>>();

        let min_bin_size = if to_refine {
            REFINING_BIN_SIZE
        } else {
            self.min_bin_size
        };

        let mut reader = parse_fastx_file(path::Path::new(&self.assembly))?;
        let mut writers: HashMap<String, BufWriter<File>> = HashMap::new();
        let mut single_contig_bin_id = 0;
        let mut unrecognised = 0;
        let mut written = 0;

        while let Some(record) = reader.next() {
            let seqrec = record?;
            let contig_name = std::str::from_utf8(seqrec.id())?.to_string();
            let contig_length = seqrec.seq().len();

            let cluster_label = if contig_length < self.min_contig_size {
                self.leftover_label(contig_length, self.min_bin_size, &mut single_contig_bin_id)
            } else {
                match labels.get(contig_name.as_str()) {
                    Some(Some(cluster_label)) => format!("{cluster_label}"),
                    Some(None) => {
                        self.leftover_label(contig_length, min_bin_size, &mut single_contig_bin_id)
                    }
                    None => {
                        unrecognised += 1;
                        self.leftover_label(contig_length, min_bin_size, &mut single_contig_bin_id)
                    }
                }
            };

            let writer = match writers.entry(cluster_label) {
                Entry::Occupied(entry) => entry.into_mut(),
                Entry::Vacant(entry) => {
                    let bin_path = path::Path::new(&self.output_directory).join(format!(
                        "rosella_bin_{}{}",
                        entry.key(),
                        RECOVER_FASTA_EXTENSION
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
        }

        if unrecognised > 0 {
            warn!(
                "{} contigs in the assembly were not in the coverage table and went unbinned",
                unrecognised
            );
        }
        debug!("Wrote {} contigs into {} bins", written, writers.len());

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

#[derive(Debug, Eq, PartialEq, PartialOrd)]
pub struct ClusterResult {
    pub(crate) contig_index: usize,
    pub(crate) cluster_label: Option<usize>,
}

impl ClusterResult {
    pub fn new(contig_index: usize, cluster_label: Option<usize>) -> Self {
        Self {
            contig_index,
            cluster_label,
        }
    }
}

impl Ord for ClusterResult {
    fn cmp(&self, other: &Self) -> Ordering {
        self.contig_index.cmp(&other.contig_index)
    }
}
