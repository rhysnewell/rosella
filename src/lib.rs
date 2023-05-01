pub mod cli;
pub mod coverage;
pub mod external;
pub mod kmers;
pub mod recover;

use log::info;

pub const AUTHOR: &str =
    "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology";
pub const AUTHOR_AND_EMAIL: &str =
    "Rhys J. P. Newell, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology <rhys.newell94 near gmail.com>";
pub const EMAIL: &str = "rhys.newell94 near gmail.com";

pub fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> f32 {
    match m.contains_id(parameter) {
        true => {
            let mut percentage: f32 = *m.get_one(parameter).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                panic!("Invalid alignment percentage: '{}'", percentage);
            }
            info!("Using {} {}%", parameter, percentage * 100.0);
            percentage
        }
        false => 0.0,
    }
}

// Enum for exclusion out here so long read can find it
pub enum GenomeExclusionTypes {
    SeparatorType,
    NoneType,
    GenomesAndContigsType,
}