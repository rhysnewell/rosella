use std::collections::HashMap;
use std::fmt::Display;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::Result;

use crate::markers::{CHECKM, COPY_BINS, ContigMarkers, GTDB, Place, Reading};
use crate::quality::Scorer;
use crate::quality::bases::{self, Bases};

pub const SETS_FILE: &str = "quality_sets.tsv";
pub const DUPLICATES_FILE: &str = "marker_duplicates.tsv";

pub struct Scored<'a> {
    pub markers: &'a ContigMarkers,
    pub names: &'a [String],
    pub lengths: &'a [usize],
    pub bases: &'a HashMap<usize, Bases>,
}

pub struct Bin<'a> {
    pub name: String,
    pub contigs: &'a [usize],
    pub strain: Option<f64>,
}

fn or_na<T: Display>(value: Option<T>) -> String {
    value.map_or_else(|| "NA".to_string(), |held| held.to_string())
}

fn sink(path: &Path) -> Result<BufWriter<std::fs::File>> {
    Ok(BufWriter::new(std::fs::File::create(path)?))
}

fn panel_header(panel: &str) -> String {
    let copies = (0..COPY_BINS)
        .map(|at| match at {
            last if last == COPY_BINS - 1 => format!("{panel}_copies_{at}plus"),
            _ => format!("{panel}_copies_{at}"),
        })
        .collect::<Vec<_>>()
        .join("\t");
    format!("{panel}_markers\t{panel}_marker_groups\t{copies}")
}

fn panel_fields(reading: Option<&Reading>) -> String {
    let Some(reading) = reading else {
        return ["NA"; COPY_BINS + 2].join("\t");
    };
    let copies = reading.copies.map(|count| count.to_string()).join("\t");
    format!("{}\t{}\t{copies}", reading.markers, reading.groups)
}

fn gene_places(genes: &[Place]) -> String {
    if genes.is_empty() {
        return "NA".to_string();
    }
    genes
        .iter()
        .map(|place| {
            let strand = if place.reverse { '-' } else { '+' };
            format!("{}-{}:{strand}", place.gene_begin, place.gene_end)
        })
        .collect::<Vec<_>>()
        .join(",")
}

impl Scored<'_> {
    pub fn write(&self, bins: &[Bin], quality: &Path) -> Result<()> {
        let directory = quality.parent().unwrap_or(Path::new("."));
        self.write_quality(bins, quality)?;
        self.write_sets(bins, &directory.join(SETS_FILE))?;
        self.write_duplicates(bins, &directory.join(DUPLICATES_FILE))
    }

    fn write_quality(&self, bins: &[Bin], path: &Path) -> Result<()> {
        let mut sink = sink(path)?;
        writeln!(
            sink,
            "bin\tcontigs\tbp\tset\tgtdb_completeness\tgtdb_contamination\t\
             checkm_completeness\tcheckm_contamination\tcheckm_strain_heterogeneity\t{}\t{}\t\
             n50\tlongest_contig\tmean_contig_length\tambiguous_bases\tgc\tgc_std\t\
             genes\tcoding_density\tmean_gene_length",
            panel_header(GTDB),
            panel_header(CHECKM),
        )?;
        for bin in bins {
            let held = self.markers.score(bin.contigs);
            let readings = self.markers.readings(bin.contigs);
            let chosen = readings
                .as_ref()
                .and_then(|(chosen, sets)| sets.get(*chosen));
            let checkm = chosen.and_then(|reading| reading.checkm.as_ref());
            let checkm_quality = match checkm {
                Some(read) => format!("{:.2}\t{:.2}", read.completeness, read.contamination),
                None => "NA\tNA".to_string(),
            };
            let bp = bin
                .contigs
                .iter()
                .map(|contig| self.lengths[*contig])
                .sum::<usize>();
            writeln!(
                sink,
                "{}\t{}\t{bp}\t{}\t{:.2}\t{:.2}\t{checkm_quality}\t{}\t{}\t{}\t{}",
                bin.name,
                bin.contigs.len(),
                self.markers.set_name(held.set),
                held.completeness,
                held.contamination,
                or_na(bin.strain.map(|strain| format!("{strain:.2}"))),
                panel_fields(chosen.map(|reading| &reading.gtdb)),
                panel_fields(checkm),
                self.bin_fields(bin.contigs),
            )?;
        }
        sink.flush()?;
        Ok(())
    }

    fn bin_fields(&self, contigs: &[usize]) -> String {
        let mut lengths = contigs
            .iter()
            .map(|contig| self.lengths[*contig])
            .collect::<Vec<_>>();
        let bp = lengths.iter().sum::<usize>();
        let longest = lengths.iter().copied().max().unwrap_or_default();
        let mean = bp.checked_div(lengths.len()).unwrap_or_default();
        let n50 = bases::n50(&mut lengths);
        let composition = bases::composition(self.bases, self.lengths, contigs);
        let coding = self.markers.coding(contigs);
        let density = coding
            .filter(|coding| coding.annotated_bp > 0)
            .map(|coding| {
                format!(
                    "{:.2}",
                    100.0 * coding.coding_bases as f64 / coding.annotated_bp as f64
                )
            });
        let gene = coding
            .filter(|coding| coding.genes > 0)
            .map(|coding| coding.coding_bases / coding.genes);
        format!(
            "{n50}\t{longest}\t{mean}\t{}\t{}\t{}\t{}\t{}\t{}",
            or_na(composition.map(|held| held.ambiguous)),
            or_na(composition.map(|held| format!("{:.2}", held.gc))),
            or_na(composition.map(|held| format!("{:.2}", held.gc_spread))),
            or_na(coding.map(|held| held.genes)),
            or_na(density),
            or_na(gene),
        )
    }

    fn write_sets(&self, bins: &[Bin], path: &Path) -> Result<()> {
        let mut sink = sink(path)?;
        writeln!(sink, "bin\tset\tpanel\tcompleteness\tcontamination\tchosen")?;
        for bin in bins {
            let Some((chosen, sets)) = self.markers.readings(bin.contigs) else {
                continue;
            };
            for held in sets {
                let name = self.markers.set_name(held.set as u16);
                let picked = u8::from(held.set == chosen);
                let panels = [(GTDB, Some(held.gtdb)), (CHECKM, held.checkm)];
                for (panel, reading) in panels {
                    let Some(reading) = reading else {
                        continue;
                    };
                    writeln!(
                        sink,
                        "{}\t{name}\t{panel}\t{:.2}\t{:.2}\t{picked}",
                        bin.name, reading.completeness, reading.contamination,
                    )?;
                }
            }
        }
        sink.flush()?;
        Ok(())
    }

    fn write_duplicates(&self, bins: &[Bin], path: &Path) -> Result<()> {
        let mut sink = sink(path)?;
        writeln!(sink, "bin\tpanel\tmarker\tcontig\tcopies\tgenes")?;
        for bin in bins {
            for row in self.markers.duplicates(bin.contigs) {
                writeln!(
                    sink,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    bin.name,
                    row.panel,
                    row.marker,
                    self.names[row.contig],
                    row.copies,
                    gene_places(&row.genes),
                )?;
            }
        }
        sink.flush()?;
        Ok(())
    }
}
