use rosella::progress::{Stage, counted, done, hold, spinning};

#[test]
fn a_held_pass_carries_its_inner_calls_and_lets_go_when_done() {
    let pass = hold(spinning(Stage::Partitioning));
    pass.set_message("pass 2");
    let inner = counted(Stage::NearestNeighbours, 30);
    inner.inc(1);
    assert_eq!(pass.message(), "pass 2, nearest neighbours 1/30");
    done(&inner);
    assert_eq!(pass.message(), "pass 2");
    done(&pass);

    let next = hold(spinning(Stage::RescuingUnbinned));
    let inner = counted(Stage::Partitioning, 10);
    inner.inc(1);
    assert_eq!(next.message(), "partitioning 1/10");
    done(&inner);
    done(&next);
}
