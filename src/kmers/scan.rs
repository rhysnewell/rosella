use anyhow::Result;
use ndarray::Array2;

use crate::kmers::halves::Halves;
use crate::kmers::kmer_counting::{KmerFrequencyTable, block_of, counts_of, frequencies_of};
use crate::kmers::sketch::{ContigSketches, SketchParams, sketch_sequence};

#[derive(Clone, Copy, Debug, Default)]
pub struct Floors {
    pub composition: Option<usize>,
    pub sketch: Option<usize>,
    pub halves: Option<usize>,
}

pub struct Scanned {
    pub composition: KmerFrequencyTable,
    pub sketches: Option<ContigSketches>,
    pub halves: Halves,
}

// Every table built from the sequences comes out of one read, because inflating the assembly
// holds a core for the whole read while the measures beside it cost little.
pub fn scan(assembly: &str, kmer_size: usize, floors: Floors) -> Result<Scanned> {
    let block = block_of(kmer_size);
    let params = SketchParams::default();
    let reaches = |floor: Option<usize>, bases: usize| floor.is_some_and(|floor| bases >= floor);
    let mut names = Vec::new();
    let mut rows = Vec::new();
    let mut sketches = floors.sketch.map(|_| ContigSketches::new());
    let mut halves = Halves::new(block.width);
    let lowest = [floors.composition, floors.sketch, floors.halves]
        .into_iter()
        .flatten()
        .min();
    if let Some(lowest) = lowest {
        let progress = crate::progress::spinning(crate::progress::Stage::CountingKmers);
        crate::kmers::measured(
            assembly,
            lowest,
            |sequence, bases| {
                let halved = reaches(floors.halves, bases).then(|| {
                    let (first, second) = sequence.split_at(sequence.len() / 2);
                    let counted = [counts_of(first, &block), counts_of(second, &block)];
                    (sequence.len(), counted)
                });
                (
                    reaches(floors.composition, bases).then(|| frequencies_of(sequence, &block)),
                    reaches(floors.sketch, bases).then(|| sketch_sequence(sequence, params)),
                    halved,
                )
            },
            |chunk| {
                for (name, measures) in chunk {
                    let (row, sketch, halved) = measures.unwrap_or((None, None, None));
                    if let Some(halved) = halved {
                        halves.push(name.clone(), halved);
                    }
                    if let Some(sketches) = &mut sketches {
                        sketches.push(name.clone(), sketch.unwrap_or_default());
                    }
                    if let Some(row) = row {
                        names.push(name);
                        rows.extend(row);
                    }
                }
                progress.set_message(format!("{} contigs", names.len()));
            },
        )?;
        progress.finish_and_clear();
    }
    Ok(Scanned {
        composition: KmerFrequencyTable::new(
            kmer_size,
            Array2::from_shape_vec((names.len(), block.width), rows)?,
            names,
        ),
        sketches,
        halves,
    })
}
