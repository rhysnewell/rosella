use std::collections::VecDeque;

use crate::clustering::codelength::plogp;
use crate::clustering::graph_partition::{Incident, compact, node_degrees, visit_order};
use crate::clustering::leiden::{Level, aggregate};
use crate::embedding::Graph;

const MAX_LEVELS: usize = 20;

fn module_term(exit: f64, visit: f64) -> f64 {
    plogp(exit + visit) - 2.0 * plogp(exit)
}

struct Modules {
    exit: Vec<f64>,
    visit: Vec<f64>,
    total_exit: f64,
    two_m: f64,
}

impl Modules {
    fn singletons(level: &Level, out: &[f64], two_m: f64) -> Self {
        Self {
            exit: out.to_vec(),
            visit: level.size.clone(),
            total_exit: out.iter().sum(),
            two_m,
        }
    }

    fn term(&self, exit: f64, visit: f64) -> f64 {
        module_term(exit / self.two_m, visit / self.two_m)
    }

    fn escape(&self, total_exit: f64) -> f64 {
        plogp(total_exit / self.two_m)
    }
}

/// Rosvall and Bergstrom (2008) applied to the moves themselves rather than to a ranking over
/// a resolution ladder, so nothing here has a resolution to set.
fn local_move(level: &Level, two_m: f64, seed: u64) -> Vec<usize> {
    let nodes = level.len();
    let out = (0..nodes)
        .map(|node| {
            level.neighbours[node]
                .iter()
                .map(|(_, weight)| *weight)
                .sum::<f64>()
        })
        .collect::<Vec<f64>>();

    let mut of = (0..nodes).collect::<Vec<_>>();
    let mut modules = Modules::singletons(level, &out, two_m);
    let mut incident = Incident::new(nodes);
    let mut queued = vec![true; nodes];
    let mut queue = visit_order(nodes, seed)
        .into_iter()
        .collect::<VecDeque<_>>();

    while let Some(node) = queue.pop_front() {
        queued[node] = false;
        let current = of[node];
        level.gather(&mut incident, node, &of);

        let weight = level.size[node];
        let held = incident.get(current);
        let exit_without = modules.exit[current] - out[node] + 2.0 * held;
        let visit_without = modules.visit[current] - weight;
        let escape_without = modules.total_exit - modules.exit[current] + exit_without;

        let stay = modules.escape(modules.total_exit)
            + modules.term(modules.exit[current], modules.visit[current]);
        let leave = modules.term(exit_without, visit_without);

        let mut best = current;
        let mut best_gain = 0.0;
        let mut best_exit = 0.0;
        for (target, shared) in incident.iter() {
            if target == current {
                continue;
            }
            let exit_with = modules.exit[target] + out[node] - 2.0 * shared;
            let escape = escape_without - modules.exit[target] + exit_with;
            let gain = stay + modules.term(modules.exit[target], modules.visit[target])
                - modules.escape(escape)
                - leave
                - modules.term(exit_with, modules.visit[target] + weight);
            if gain > best_gain || (gain == best_gain && gain > 0.0 && target < best) {
                best = target;
                best_gain = gain;
                best_exit = exit_with;
            }
        }

        if best == current {
            continue;
        }

        modules.total_exit = escape_without - modules.exit[best] + best_exit;
        modules.exit[current] = exit_without;
        modules.visit[current] = visit_without;
        modules.exit[best] = best_exit;
        modules.visit[best] += weight;
        of[node] = best;

        for (target, _) in &level.neighbours[node] {
            if of[*target] != best && !queued[*target] {
                queued[*target] = true;
                queue.push_back(*target);
            }
        }
    }

    of
}

/// The map equation is what ranks a labelling, so this is the search that optimises the score
/// it is judged on. Needs no resolution, and so has no ladder to select a rung from.
pub fn infomap(graph: &Graph, seed: u64) -> Vec<i32> {
    let nodes = graph.rows();
    let degrees = node_degrees(graph);
    let two_m = degrees.iter().sum::<f64>();
    if nodes == 0 || two_m <= 0.0 {
        return vec![0; nodes];
    }

    let mut level = Level::from_graph(graph).with_size(degrees);
    let mut membership = (0..nodes).collect::<Vec<_>>();
    let mut labels = vec![0i32; nodes];

    for _ in 0..MAX_LEVELS {
        let of = local_move(&level, two_m, seed);
        for (node, slot) in labels.iter_mut().enumerate() {
            *slot = of[membership[node]] as i32;
        }

        let (next, ids) = aggregate(&level, &of);
        if next.len() == level.len() {
            break;
        }
        for slot in &mut membership {
            *slot = ids[*slot];
        }
        level = next;
    }

    compact(&labels)
}
