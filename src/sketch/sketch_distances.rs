use std::cmp::Ordering;
use hnsw_rs::dist::Distance;
use finch::{serialization::Sketch, sketch_schemes::KmerCount};
use serde_derive::{Serialize, Deserialize};

pub fn distance(
    query_sketch: &Sketch,
    ref_sketch: &Sketch,
) -> Result<SketchDistanceContainer, &'static str> {
    // since we always examine to the lowest of the sketch maxima, a
    // min_scale of 0 is a noop; otherwise we only set a scale if both of
    // the sketches are scaled (there may be a slight improvement in
    // comparing a unscaled "higher range" sketch to a scaled lower range
    // using the scale, but that makes things more complicated because we
    // need two scale values, etc)
    let mut min_scale = 0.;
    if let Some(scale1) = query_sketch.sketch_params.hash_info().3 {
        if let Some(scale2) = ref_sketch.sketch_params.hash_info().3 {
            min_scale = f64::min(scale1, scale2);
        }
    }
    let distances = raw_distance(&query_sketch.hashes, &ref_sketch.hashes, min_scale);


    let containment = distances.0;
    let jaccard = distances.1;
    let min_jaccard = distances.2;
    let common_hashes = distances.3;
    let total_hashes = distances.4;
    let k = query_sketch.sketch_params.k() as f64;
    let mash_distance: f64 = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / k;
    Ok(SketchDistanceContainer {
        containment,
        jaccard,
        min_jaccard,
        mash_distance: f64::min(1f64, f64::max(0f64, mash_distance)),
        common_hashes,
        total_hashes,
        query: query_sketch.name.to_string(),
        reference: ref_sketch.name.to_string(),
    })
}

/// Estimates set statistics based on two slices of `KmerCount` sorted by hash,
/// ignoring hashes past a certain point (details below).
///
/// Returns a tuple of the form (containment,
/// [jaccard index](https://en.wikipedia.org/wiki/Jaccard_index),
/// size of [intersection](https://en.wikipedia.org/wiki/Intersection_(set_theory)),
/// size of [union](https://en.wikipedia.org/wiki/Union_(set_theory))).
///
/// Ignores hashes that are greater than `scale` * the maximum possible hash
/// (i.e. `u64::max_value()`). A scale of `0_f64` indicates to not ignore
/// hashes for this reason. Additionally ignores a slice's hashes if they are
/// greater than the other slice's maximum hash; this method better
/// approximates the document distance when the documents are of different
/// or of unknown sizes.
///
/// If the `KmerCount` slice arguments are not sorted by hash, the values of
/// the returned set statistics are unspecified.
pub fn raw_distance(
    query_hashes: &[KmerCount],
    ref_hashes: &[KmerCount],
    scale: f64,
) -> (f64, f64, f64, u64, u64) {

    let mut i: usize = 0;
    let mut j: usize = 0;
    let mut common: u64 = 0;
    while let (Some(query), Some(refer)) = (query_hashes.get(i), ref_hashes.get(j)) {
        match query.hash.cmp(&refer.hash) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                common += 1;
                i += 1;
                j += 1;
            }
        }
    }

    // at this point we've exhausted one of the two sketches, but we may have
    // more counts in the other to compare if these were scaled sketches
    if scale > 0. {
        let max_hash = u64::max_value() / scale.recip() as u64;
        while query_hashes
            .get(i)
            .map(|kmer_count| kmer_count.hash < max_hash)
            .unwrap_or(false)
        {
            i += 1;
        }
        while ref_hashes
            .get(j)
            .map(|kmer_count| kmer_count.hash < max_hash)
            .unwrap_or(false)
        {
            j += 1;
        }
    }

    let containment = if j == 0 { 0. } else { common as f64 / j as f64 };
    let total = i as u64 - common + j as u64;
    let jaccard: f64 = if total == 0 {
        1.
    } else {
        common as f64 / total as f64
    };

    let min_jaccard: f64 = if i == 0 && j == 0 {
        1.
    } else {
        common as f64 / std::cmp::min(i, j) as f64
    };

    (containment, jaccard, min_jaccard, common, total)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SketchDistanceContainer {
    pub containment: f64,
    pub jaccard: f64,
    pub min_jaccard: f64,
    #[serde(rename = "mashDistance")]
    pub mash_distance: f64,
    #[serde(rename = "commonHashes")]
    pub common_hashes: u64,
    #[serde(rename = "totalHashes")]
    pub total_hashes: u64,
    pub query: String,
    pub reference: String,
}

pub struct SketchDistance;
impl Distance<&Sketch> for SketchDistance {
    fn eval(&self, va: &[&Sketch], vb: &[&Sketch]) -> f32 {
        // take the mean of the distances between all pairs of sketches
        let mut result = (va.iter()
            .map(|query_sketch| {
                vb.iter()
                    .map(|ref_sketch| 1.0 - distance(query_sketch, ref_sketch).unwrap().min_jaccard)
                    .sum::<f64>()
                    / vb.len() as f64
            })
            .sum::<f64>()
            / va.len() as f64)
            as f32;
        
        if result.is_nan() {
            result = 1.0;
        }
        result
    }
}