use crate::markers::MarkerSet;
use crate::markers::hmm_table::{Best, Reach};
use crate::quality::orfs::Orf;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Hit {
    pub marker: u16,
    pub partial: bool,
    pub place: Place,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Place {
    pub reach: Reach,
    pub score: f32,
    pub gene_begin: u32,
    pub gene_end: u32,
    pub reverse: bool,
    pub cut_left: bool,
    pub cut_right: bool,
}

const FIELDS: usize = 12;

impl Hit {
    pub fn called(marker: u16, orf: &Orf, best: &Best) -> Self {
        Self {
            marker,
            partial: orf.partial(),
            place: Place {
                reach: best.reach,
                score: best.score as f32,
                gene_begin: orf.begin,
                gene_end: orf.begin + orf.bases as u32 - 1,
                reverse: orf.reverse,
                cut_left: orf.cut_left,
                cut_right: orf.cut_right,
            },
        }
    }

    // Two halves of one gene cut by a contig end align to different ends of the model, while a
    // second copy aligns over the same stretch.
    pub fn same_part(&self, other: &Self) -> bool {
        let (ours, theirs) = (self.place.reach, other.place.reach);
        let covered = |reach: Reach| (reach.model_to + 1).saturating_sub(reach.model_from);
        let shared = (ours.model_to.min(theirs.model_to) + 1)
            .saturating_sub(ours.model_from.max(theirs.model_from));
        2 * shared > covered(ours).min(covered(theirs))
    }

    pub fn fields(&self, separator: char) -> String {
        let place = &self.place;
        let reach = place.reach;
        [
            reach.model_from.to_string(),
            reach.model_to.to_string(),
            reach.model_length.to_string(),
            reach.protein_from.to_string(),
            reach.protein_to.to_string(),
            place.score.to_string(),
            place.gene_begin.to_string(),
            place.gene_end.to_string(),
            if place.reverse { "-" } else { "+" }.to_string(),
            u8::from(place.cut_left).to_string(),
            u8::from(place.cut_right).to_string(),
        ]
        .join(&separator.to_string())
    }

    pub fn encode(&self, set: &MarkerSet) -> String {
        format!("{}:{}", set.name(self.marker), self.fields(':'))
    }

    pub fn decode(field: &str, set: &MarkerSet) -> Option<Self> {
        let parts = field.split(':').collect::<Vec<_>>();
        if parts.len() != FIELDS {
            return None;
        }
        let number = |at: usize| parts[at].parse::<u32>().ok();
        let flag = |at: usize| match parts[at] {
            "1" => Some(true),
            "0" => Some(false),
            _ => None,
        };
        let (cut_left, cut_right) = (flag(10)?, flag(11)?);
        Some(Self {
            marker: set.id(parts[0])?,
            partial: cut_left || cut_right,
            place: Place {
                reach: Reach {
                    model_from: number(1)?,
                    model_to: number(2)?,
                    model_length: number(3)?,
                    protein_from: number(4)?,
                    protein_to: number(5)?,
                },
                score: parts[6].parse().ok()?,
                gene_begin: number(7)?,
                gene_end: number(8)?,
                reverse: parts[9] == "-",
                cut_left,
                cut_right,
            },
        })
    }
}
