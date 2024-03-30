use crate::pathwise_graph::PathGraph;

#[derive(Debug, Clone)]
pub struct DisplacementMatrix {
    pub dfs: Vec<isize>,
    pub dfe: Vec<isize>,
}

impl DisplacementMatrix {
    pub fn new(graph: &PathGraph, rev_graph: &PathGraph) -> Self {
        DisplacementMatrix {
            dfs: get_distance_from_start(rev_graph),
            dfe: get_distance_from_end(graph),
        }
    }

    pub fn get_displ(&self, node_i: usize, node_j: usize) -> i32 {
        return if node_i == node_j {
            0
        } else {
            ((self.dfs[node_i] - self.dfs[node_j]).abs()
                + (self.dfe[node_i] - self.dfe[node_j]).abs()) as i32
        };
    }
}

fn get_distance_from_start(graph: &PathGraph) -> Vec<isize> {
    let nwp = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let lnz_len = graph.lnz.len();
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[0] = 0;
    for (p, _) in pred_hash.get_preds_and_paths(0) {
        r_values[p] = 1;
    }

    for i in 1..lnz_len - 1 {
        if r_values[i] == -1 || r_values[i] > r_values[i - 1] + 1 {
            r_values[i] = r_values[i - 1] + 1;
        }
        if nwp[i] {
            for (p, _) in pred_hash.get_preds_and_paths(i) {
                if r_values[p] == -1 || r_values[p] > r_values[i] + 1 {
                    r_values[p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values
}
fn get_distance_from_end(graph: &PathGraph) -> Vec<isize> {
    let nwp = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let lnz_len = graph.lnz.len();
    let mut r_values: Vec<isize> = vec![-1; lnz_len];
    r_values[lnz_len - 1] = 0;

    for (p, _) in pred_hash.get_preds_and_paths(lnz_len - 1) {
        r_values[p] = 1;
    }

    for i in (1..lnz_len - 1).rev() {
        if r_values[i] == -1 || r_values[i] > r_values[i + 1] + 1 {
            r_values[i] = r_values[i + 1] + 1;
        }
        if nwp[i] {
            for (p, _) in pred_hash.get_preds_and_paths(i) {
                if r_values[p] == -1 || r_values[p] > r_values[i] + 1 {
                    r_values[p] = r_values[i] + 1;
                }
            }
        }
    }
    r_values
}
