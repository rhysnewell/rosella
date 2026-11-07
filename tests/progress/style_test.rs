// `progress::styled` falls back to the default bar on a bad template, so a typo in one would
// silently lose that stage's colour rather than fail.

use indicatif::ProgressStyle;
use rosella::progress::{Stage, counted_template, spinning_template, stepped_template};

#[test]
fn every_stage_template_parses() {
    for stage in Stage::ALL {
        for template in [
            counted_template(stage),
            spinning_template(stage),
            stepped_template(stage),
        ] {
            assert!(
                ProgressStyle::with_template(&template).is_ok(),
                "{stage:?}: {template}"
            );
        }
    }
}
