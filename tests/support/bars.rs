#![allow(dead_code)]

use rosella::quality::Worth;
use rosella::refine::rung::Bars;

/// The completeness the engine hands the pool is already through the scorer's own offset, so a
/// test that wants the shipped ladder passes 80 rather than the 90 a user asks for.
pub fn bars(completeness: f64) -> Bars {
    Bars {
        min_bin_size: 0,
        completeness,
        contamination: 5.0,
        worth: Worth {
            contamination: 2.0,
            allowance: 0.0,
        },
        rung_floor: 0.56,
        rung_ceiling: 3.0,
    }
}
