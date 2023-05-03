use anyhow::Result;
use finch::{serialization::Sketch, sketch_schemes::SketchParams, filtering::FilterParams};
use ndarray::Array2;
use needletail::Sequence;
use rayon::prelude::*;

use super::sketch_distances::distance;

pub fn sketch_contigs(m: &clap::ArgMatches) -> Result<ContigSketchResult> {
    let mut sketcher = ContigSketcher::new(m);
    sketcher.run()
}

struct ContigSketcher {
    assembly: String,
    sketch_settings: SketchSettings,
}

impl ContigSketcher {
    pub fn new(m: &clap::ArgMatches) -> Self {
        let assembly = m.get_one::<String>("assembly").unwrap().clone();
        let sketch_settings = SketchSettings::new(m);

        Self {
            assembly,
            sketch_settings,
        }
    }

    pub fn run(&mut self) -> Result<ContigSketchResult> {
        let (contig_names, contig_sketches) = self.obtain_sketches()?;
        let contig_sketch_similarity = self.calculate_similarity(contig_sketches)?;

        let contig_sketch_result = ContigSketchResult {
            contig_names,
            contig_sketch_similarity,
        };

        Ok(contig_sketch_result)
    }

    fn obtain_sketches(&self) -> Result<(Vec<String>, Vec<Sketch>)> {
        let mut reader = needletail::parse_fastx_file(&self.assembly)?;
        let mut contig_names = Vec::new();
        let mut contig_sketches = Vec::new();
        while let Some(record) = reader.next() {
            let seqrec = record?;
            let contig_name = std::str::from_utf8(seqrec.id())?.to_string();
            // normalise sequence
            let seqrec = seqrec.normalize(false);
            contig_names.push(contig_name.clone());
            // get genome size
            let contig_size = seqrec.len();


            let kmers_to_sketch = (contig_size as f64 * self.sketch_settings.sketch_scale) as usize;
            let sketch_params = SketchParams::Scaled {
                kmers_to_sketch,
                kmer_length: self.sketch_settings.kmer_size as u8,
                scale: self.sketch_settings.sketch_scale,
                hash_seed: self.sketch_settings.seed,
            };

            let mut sketcher = sketch_params.create_sketcher();
            sketcher.process(&seqrec);

            let (seq_length, num_valid_kmers) = sketcher.total_bases_and_kmers();
            let hashes = sketcher.to_vec();

            // do filtering
            let sketch = Sketch {
                name: contig_name,
                seq_length,
                num_valid_kmers,
                comment: "".to_string(),
                hashes,
                filter_params: FilterParams::default(),
                sketch_params,
            };
            contig_sketches.push(sketch);
        }
        Ok((contig_names, contig_sketches))
    }

    fn calculate_similarity(&self, signatures: Vec<Sketch>) -> Result<Array2<f32>> {
        // calculate number of pairwise combinations of signatures
        let n_signatures = signatures.len();
        // efficiently calculate pairwise distances in self.signatures
        // distances is a condensed pairwise matrix
        // let mut distances = vec![1.0; n_combinations];
        let distances = (0..n_signatures)
            .into_par_iter()
            .flat_map(|i| {
                let sig_i = &signatures[i];
                ((i + 1)..n_signatures)
                    .into_par_iter()
                    .map(|j| {
                        let sig_j = &signatures[j];
                        let dist = distance(sig_i, sig_j)
                            .expect("Failed to generate distances");
                        dist.jaccard as f32
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // convert distances to a square matrix
        let mut distance_matrix = vec![vec![1.0; n_signatures]; n_signatures];

        let mut index = 0;
        for i in 0..n_signatures {
            for j in (i + 1)..n_signatures {
                distance_matrix[i][j] = distances[index];
                distance_matrix[j][i] = distances[index];
                index += 1;
            }
        }

        let distance_matrix = Array2::from_shape_vec((n_signatures, n_signatures), distances)?;
        Ok(distance_matrix)
    }
}

struct SketchSettings {
    kmer_size: usize,
    sketch_scale: f64,
    seed: u64
}

impl SketchSettings {
    pub fn new(m: &clap::ArgMatches) -> Self {
        let kmer_size = *m.get_one::<usize>("kmer-size").unwrap();
        let sketch_scale = *m.get_one::<f64>("sketch-scale").unwrap();
        let seed = *m.get_one::<u64>("seed").unwrap();
        Self {
            kmer_size,
            sketch_scale,
            seed,
        }
    }
}

pub struct ContigSketchResult {
    pub contig_names: Vec<String>,
    pub contig_sketch_similarity: Array2<f32>,
}