use crate::clustering::graph_partition::visit_order;
use crate::clustering::sbm::state::BlockState;

const MOVE_ROUNDS: usize = 10;
const MERGE_RATIO: f64 = 2.0;

fn gather(out: &mut Vec<(usize, f64)>, edges: &[(usize, f64)], of: &[usize]) {
    out.clear();
    out.extend(edges.iter().map(|(target, weight)| (of[*target], *weight)));
    out.sort_unstable_by_key(|(block, _)| *block);
    out.dedup_by(|(later, weight), (kept, total)| {
        if later == kept {
            *total += *weight;
            true
        } else {
            false
        }
    });
}

fn moves(state: &mut BlockState, neighbours: &[Vec<(usize, f64)>], node_degree: &[f64], seed: u64) {
    let mut incidence = Vec::new();
    for round in 0..MOVE_ROUNDS {
        let mut moved = false;
        for node in visit_order(neighbours.len(), seed.wrapping_add(round as u64)) {
            gather(&mut incidence, &neighbours[node], &state.of);
            if incidence.is_empty() {
                continue;
            }
            let from = state.of[node];
            let mut best = from;
            let mut best_cost = 0.0;
            for (block, _) in &incidence {
                if *block == from {
                    continue;
                }
                let cost = state.move_cost(node, *block, &incidence, node_degree[node]);
                if cost < best_cost {
                    best = *block;
                    best_cost = cost;
                }
            }
            if best != from {
                state.apply_move(node, best, &incidence, node_degree[node]);
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }
}

/// Costs are read off the state before any of the round's merges land, so a block that has just
/// taken one is held back rather than merged again on a stale number.
fn merge_down(state: &mut BlockState, target: usize) -> bool {
    let mut candidates = state
        .live_blocks()
        .into_iter()
        .filter_map(|gone| {
            state
                .neighbouring_blocks(gone)
                .filter(|keep| *keep != gone && state.count(*keep) > 0)
                .map(|keep| (state.merge_cost(keep, gone), keep, gone))
                .min_by(|a, b| a.0.total_cmp(&b.0))
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut redirect = (0..state.capacity()).collect::<Vec<_>>();
    let mut held = vec![false; state.capacity()];
    let mut merged = 0;

    for (_, keep, gone) in candidates {
        if state.blocks() <= target {
            break;
        }
        if held[gone] || held[keep] || state.count(gone) == 0 || state.count(keep) == 0 {
            continue;
        }
        state.apply_merge(keep, gone);
        redirect[gone] = keep;
        held[gone] = true;
        held[keep] = true;
        merged += 1;
    }

    state.relabel(&redirect);
    merged > 0
}

/// Node moves against the description length, then agglomeration down through the block counts,
/// keeping whichever count described the graph most cheaply (Peixoto, 2014).
pub fn search(
    neighbours: &[Vec<(usize, f64)>],
    node_degree: &[f64],
    of: Vec<usize>,
    block_count: usize,
    seed: u64,
) -> Vec<usize> {
    let mut state = BlockState::new(neighbours, node_degree, of, block_count);
    let mut best = state.of.clone();
    let mut best_length = f64::INFINITY;

    loop {
        moves(&mut state, neighbours, node_degree, seed);
        let length = state.description_length();
        if length < best_length {
            best_length = length;
            best = state.of.clone();
        }
        if state.blocks() <= 1 {
            break;
        }
        let target = ((state.blocks() as f64 / MERGE_RATIO) as usize).max(1);
        if !merge_down(&mut state, target) {
            break;
        }
    }

    best
}
