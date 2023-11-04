use std::{collections::{HashSet, HashMap}, process::Command, io::{BufRead, Read, BufWriter}, cmp::max, path, fs::OpenOptions};

use anyhow::Result;
use bird_tool_utils::clap_utils::parse_list_of_genome_fasta_files;
use log::{info, debug, error};
use needletail::{parse_fastx_file, parser::{LineEnding, write_fasta}};
use rayon::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};

use crate::{coverage::{coverage_table::CoverageTable, coverage_calculator::calculate_coverage}, external_command_checker::check_for_flight, recover::recover_engine::{ClusterResult, RECOVER_FASTA_EXTENSION, UNBINNED}, kmers::kmer_counting::count_kmers};

pub const UNCHANGED_LOC: &str = "unchanged_bins";
pub const UNCHANGED_BIN_TAG: &str = "unchanged";
pub const REFINED_LOC: &str = "refined_bins";

pub fn run_refine(m: &clap::ArgMatches) -> Result<()> {
    let mut refine_engine = RefineEngine::new(m)?;
    refine_engine.run(m.get_one::<String>("bin-tag").unwrap().as_str())
}

pub struct RefineEngine {
    pub(crate) assembly: String,
    pub(crate) output_directory: String,
    pub(crate) threads: usize,
    pub(crate) kmer_frequencies: String,
    pub(crate) coverage_table: CoverageTable,
    pub(crate) checkm_results: Option<String>,
    pub(crate) mags_to_refine: Vec<String>,
    pub(crate) min_contig_size: usize,
    pub(crate) min_bin_size: usize,
    pub(crate) n_neighbours: usize,
    pub(crate) bin_unbinned: bool,
}

impl RefineEngine {
    fn new(m: &clap::ArgMatches) -> Result<Self> {
        let assembly = m.get_one::<String>("assembly").unwrap().clone();
        let output_directory = m.get_one::<String>("output-directory").unwrap().clone();
        let threads = *m.get_one::<usize>("threads").unwrap();
        let coverage_table = calculate_coverage(m)?;
        let n_contigs = coverage_table.contig_names.len();

        let kmer_frequencies = if let Some(kmer_table_path) = m.get_one::<String>("kmer-frequency-file") {
            info!("Reading TNF table.");
            kmer_table_path.clone()
        } else {
            info!("Calculating TNF table.");
            let kmer_table = count_kmers(m, Some(n_contigs))?;
            let kmer_table_path = kmer_table.table_path;
            kmer_table_path
        };
        
        let checkm_results = match m.get_one::<String>("checkm-results") {
            Some(checkm_results) => Some(checkm_results.clone()),
            None => None,
        
        };
        let mags_to_refine = match parse_list_of_genome_fasta_files(m, true) {
            Ok(mags) => mags,
            Err(e) => {
                bail!("Failed to parse MAGs to refine: {}", e);
            }
        };

        let min_contig_size = *m.get_one::<usize>("min-contig-size").unwrap();
        let min_bin_size = *m.get_one::<usize>("min-bin-size").unwrap();
        let n_neighbours = *m.get_one::<usize>("n-neighbours").unwrap();


        Ok(Self {
            assembly,
            output_directory,
            threads,
            kmer_frequencies,
            coverage_table,
            checkm_results,
            mags_to_refine,
            min_contig_size,
            min_bin_size,
            n_neighbours,
            bin_unbinned: false,
        })
    }


    pub fn run(&mut self, refined_bin_tag: &str) -> Result<()> {
        check_for_flight()?;
        // ensure output directory exists
        std::fs::create_dir_all(&self.output_directory)?;
        std::fs::create_dir_all(format!("{}/{}", &self.output_directory, UNCHANGED_LOC))?;
        std::fs::create_dir_all(format!("{}/{}", &self.output_directory, REFINED_LOC))?;
        
        let removal_string = if self.bin_unbinned {
            "small_unbinned"
        } else {
            "unbinned"
        };

        // check if we have unbinned in self.mags_to_refine
        // if so move them straigh to unchanged
        self.mags_to_refine.iter().for_each(|genome| {
            if genome.contains(removal_string) {
                self.copy_bin_to_output(genome, UNCHANGED_LOC, UNCHANGED_BIN_TAG).unwrap();
            }
        });

        self.mags_to_refine.retain(|genome| !genome.contains(removal_string));
        

        let extra_threads = max(self.threads / self.mags_to_refine.len(), 1);

        info!("Beginning refinement of {} MAGs", self.mags_to_refine.len());
        let progress_bar = ProgressBar::new(self.mags_to_refine.len() as u64);
        progress_bar.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}"));
            // .progress_chars("##-"));
        progress_bar.set_message("Refining MAGs");

        let mut results = self.mags_to_refine
            .par_iter()
            .filter(|genome| !genome.contains(UNBINNED))
            .map(|genome| {
            let result = self.run_flight_refine(genome, extra_threads);
            progress_bar.inc(1);
            (result, genome)
        }).collect::<Vec<(_, _)>>();

        if self.bin_unbinned {
            // bin uninbinned with all threads
            let unbinned_results = self.mags_to_refine
                .iter()
                .filter(|genome| genome.contains(UNBINNED))
                .map(|genome| {
                    let result = self.run_flight_refine(genome, self.threads);
                    progress_bar.inc(1);
                    (result, genome)
                }).collect::<Vec<(_, _)>>();
            
            results.extend(unbinned_results);
        }
        progress_bar.finish_with_message("Finished refining MAGs");

        // read in the results
        let mut unified_cluster_map: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (result, genome_path) in results {
            let result = result?;
            debug!("Reading in results from {} for genome {}", &result, &genome_path);
            let cluster_map = self.read_in_results(&result, genome_path)?;
            for (cluster, contigs) in cluster_map {
                // check if the cluster already exists, if so assign new cluster id
                let cluster = match unified_cluster_map.get(&cluster) {
                    Some(_) => {
                        // max key value
                        let mut max_cluster = *unified_cluster_map.keys().max().unwrap();
                        max_cluster += 1;
                        max_cluster
                    },
                    None => cluster
                };

                unified_cluster_map.insert(cluster, contigs);
            }
        }

        // get outliers
        let outliers = unified_cluster_map.remove(&0);
        let cluster_results = self.get_cluster_result(unified_cluster_map, outliers.unwrap_or_else(|| HashSet::new()));

        // write the bins
        self.write_clusters(REFINED_LOC, cluster_results, refined_bin_tag)?;
        // remove all excess JSON files in output_directory
        let excess_json_files = std::fs::read_dir(&self.output_directory)?
            .filter_map(|entry| {
                let entry = entry.unwrap();
                let path = entry.path();
                if path.extension().unwrap_or_default() == "json" {
                    Some(path)
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        for path in excess_json_files {
            std::fs::remove_file(path)?;
        }

        Ok(())
    }

    fn run_flight_refine(&self, genome_path: &str, extra_threads: usize) -> Result<String> {

        // get output prefix from genome path by remove path and extensions
        let output_prefix = genome_path
            .split("/")
            .collect::<Vec<_>>()
            .last().unwrap().to_string();

        let mut flight_cmd = Command::new("flight");
        flight_cmd.arg("refine");
        // flight_cmd.arg("--assembly").arg(&self.assembly);
        flight_cmd.arg("--genome_paths").arg(genome_path);
        flight_cmd.arg("--input").arg(&self.coverage_table.output_path);
        flight_cmd.arg("--kmer_frequencies").arg(&self.kmer_frequencies);
        flight_cmd.arg("--output_directory").arg(&self.output_directory);
        flight_cmd.arg("--output_prefix").arg(&output_prefix);
        flight_cmd.arg("--min_contig_size").arg(format!("{}", self.min_contig_size));
        flight_cmd.arg("--min_bin_size").arg(format!("{}", self.min_bin_size));
        flight_cmd.arg("--n_neighbors").arg(format!("{}", self.n_neighbours));
        flight_cmd.arg("--cores").arg(format!("{}", extra_threads));
        if let Some(checkm_results) = &self.checkm_results {
            flight_cmd.arg("--checkm_file").arg(checkm_results);
            flight_cmd.arg("--contaminated_only");
        }


        flight_cmd.stdout(std::process::Stdio::piped());
        flight_cmd.stderr(std::process::Stdio::piped());

        let mut child = match flight_cmd.spawn() {
            Ok(child) => child,
            Err(e) => {
                bail!("Error running flight: {}", e);
            }
        };

        // get the exit code
        let exit_status = child.wait()?;
        if !exit_status.success() {
            if let Some(stderr) = child.stderr.take() {
                let stderr = std::io::BufReader::new(stderr);
                for line in stderr.lines() {
                    let line = line?;
                    let message = line.split("INFO: ").collect::<Vec<_>>();
                    error!("{}", message[message.len() - 1]);
                }
            }
            bail!("Flight failed with exit code: {}", exit_status);
        }
        
        let output_json = format!("{}/{}.json", self.output_directory, output_prefix);

        Ok(output_json)
    }

    fn get_original_contig_count(&self, bin_path: &str) -> Result<usize> {
        let mut reader = parse_fastx_file(path::Path::new(&bin_path))?;
        let mut contig_count = 0;
        while let Some(_) = reader.next() {   
            contig_count += 1;
        }

        Ok(contig_count)
    }

    fn read_in_results(&self, json_path: &str, original_bin_path: &str) -> Result<HashMap<usize, HashSet<usize>>> {
        let mut bins =
            std::fs::File::open(&json_path)?;
        let mut data = String::new();
        bins.read_to_string(&mut data).unwrap();
        let mut cluster_map: HashMap<usize, HashSet<usize>> = serde_json::from_str(&data).unwrap();

        let og_contig_count = self.get_original_contig_count(original_bin_path)?;
        let mut bin_unchanged = false;
        if cluster_map.keys().len() == 1 {
            // check first values and see if length is the same as the original bin contig count
            let first_value = cluster_map.keys().next().unwrap();
            if cluster_map.get(first_value).unwrap().len() == og_contig_count {
                bin_unchanged = true;
            }
        }

        if cluster_map.len() == 0 || bin_unchanged {
            // bin was not changed, so we keep the original
            self.copy_bin_to_output(original_bin_path, UNCHANGED_LOC, UNCHANGED_BIN_TAG)?;
            
            return Ok(HashMap::new());
        }
        debug!("Cluster map: {:?} genome {}", &cluster_map, original_bin_path);
        let mut removed_bins = Vec::new();
        for (bin, contigs) in cluster_map.iter() {
            let mut contigs = contigs.iter().cloned().collect::<Vec<_>>();
            contigs.sort_unstable();

            let bin_size = contigs.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
            if bin_size < self.min_bin_size {
                removed_bins.push(*bin);
            }
        }

        for bin in removed_bins {
            match cluster_map.remove(&bin) {
                Some(contigs) => {
                    // add them to cluster 0
                    cluster_map.entry(0).or_insert_with(HashSet::new).extend(contigs.into_iter());
                },
                None => {
                    bail!("Bin {} does not exist in cluster map", bin);
                }
            }
        }

        Ok(cluster_map)
    }

    fn copy_bin_to_output(&self, genome_path: &str, bin_dir: &str, _bin_tag: &str) -> Result<()> {
        // copy to unchanged_bins
        let og_bin_name = genome_path.split("/").collect::<Vec<_>>().last().unwrap().to_string();
        // let og_bin_name_no_ext = og_bin_name.split(".").collect::<Vec<_>>()[0];
        // let new_bin_name = format!("{}.{}", og_bin_name_no_ext, RECOVER_FASTA_EXTENSION);

        let output_path = format!("{}/{}/{}", &self.output_directory, bin_dir, og_bin_name);
        std::fs::copy(&genome_path, &output_path)?;
        Ok(())
    }

    fn get_cluster_result(
        &self,
        cluster_map: HashMap<usize, HashSet<usize>>,
        outliers: HashSet<usize>
    ) -> Vec<ClusterResult> {

        let mut cluster_results = Vec::with_capacity(self.coverage_table.contig_names.len());
        for (cluster_label, contig_indices) in cluster_map.into_iter() {
            
            // check the size of the cluster and if it is too small, set the cluster label to None
            let bin_size = contig_indices.iter().map(|i| self.coverage_table.contig_lengths[*i]).sum::<usize>();
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

    fn write_clusters(&self, bin_dir: &str, cluster_results: Vec<ClusterResult>, bin_tag: &str) -> Result<()> {
        // open the assembly
        let mut reader = parse_fastx_file(path::Path::new(&self.assembly))?;

        let mut contig_idx = 0;
        let mut single_contig_bin_id = 0;
        let mut total_contig_idx = 0;
        debug!("Cluster result length: {}", cluster_results.len());
        while let Some(record) = reader.next() {
            let seqrec = record?;
            
            let contig_name = seqrec.id();
            let contig_length = seqrec.seq().len();
            if contig_length < self.min_contig_size {
                total_contig_idx += 1;
                continue;
            }


            // find cluster result that matches contig index
            let cluster_result = match cluster_results.binary_search_by_key(&total_contig_idx, |cluster_result| cluster_result.contig_index) {
                Ok(cluster_result) => &cluster_results[cluster_result],
                Err(_) => {
                    total_contig_idx += 1;
                    continue;
                }
            };

            if total_contig_idx != cluster_result.contig_index {
                debug!("Contig index mismatch. {} != {}", total_contig_idx, cluster_result.contig_index);
            }
            
            // open writer for bin
            let cluster_label = match cluster_result.cluster_label {
                Some(cluster_label) => format!("{}", cluster_label),
                None => {
                    // contig is an outlier, so check it's length. If it is greater than
                    // the minimum bin size, then bin it, otherwise discard it.
                    if seqrec.seq().len() < self.min_bin_size {
                        UNBINNED.to_string()
                    } else {
                        single_contig_bin_id += 1;
                        format!("single_contig_{}_{}", bin_tag, single_contig_bin_id)
                    }

                }
            };

            let bin_path = path::Path::new(&self.output_directory)
                .join(bin_dir)
                .join(format!("rosella_{}_{}{}", bin_tag, cluster_label, RECOVER_FASTA_EXTENSION));
            let file = OpenOptions::new().append(true).create(true).open(bin_path)?;
            let mut writer = BufWriter::new(file);
            // write contig to bin
            write_fasta(contig_name, &seqrec.seq(), &mut writer, LineEnding::Unix)?;
            contig_idx += 1;
            total_contig_idx += 1;
        }

        debug!("Final contig idx: {}", contig_idx);

        Ok(())
    }
}