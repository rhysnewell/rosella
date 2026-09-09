//! Gene calling, which the gene family scorer reads before it can score anything.

use rosella::quality::orfs::translate;

#[test]
fn translation_rewrites_the_start_and_drops_the_stop() {
    assert_eq!(translate(b"GTGAAATAA", true), "MK");
    assert_eq!(translate(b"GTGAAATAA", false), "VK");
    assert_eq!(translate(b"ATGNNNTAA", true), "MX");
}
