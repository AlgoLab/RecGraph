use std::collections::HashMap;

use crate::pathwise_graph::PathGraph;

pub fn exec(sequence: &[char], graph: &PathGraph, score_matrix: &HashMap<(char, char), i32>) {
    let lnz = &graph.lnz;
    let nodes_with_pred = &graph.nwp;
    let pred_hash = &graph.pred_hash;
    let paths_nodes = &graph.paths_nodes;
    let paths_number = graph.paths_number;
    let mut alphas = vec![paths_number + 1; lnz.len()];
    let mut dpm = vec![vec![vec![0; paths_number]; sequence.len()]; lnz.len()];

    alphas[0] = 0;
    for i in 1..lnz.len() - 1 {
        if !nodes_with_pred[i] {
            //single predecessor
            let mut common_paths = paths_nodes[i].clone();
            common_paths.and(&paths_nodes[i - 1]);

            if common_paths[alphas[i - 1]] {
                alphas[i] = alphas[i - 1];
                common_paths.iter().enumerate().for_each(|(path, is_in)| {
                    if is_in {
                        if path == alphas[i] {
                            dpm[i][0][path] =
                                dpm[i - 1][0][path] + score_matrix.get(&(lnz[i], '-')).unwrap();
                        } else {
                            dpm[i][0][path] = dpm[i - 1][0][path];
                        }
                    }
                });
            } else {
                let mut first = true;
                common_paths.iter().enumerate().for_each(|(path, is_in)| {
                    if is_in {
                        if first {
                            first = false;
                            alphas[i] = path;
                            dpm[i][0][alphas[i]] = dpm[i - 1][0][alphas[i]]
                                + dpm[i - 1][0][alphas[i - 1]]
                                + score_matrix.get(&(lnz[i], '-')).unwrap();
                        } else {
                            dpm[i][0][path] = dpm[i - 1][0][path] - dpm[i - 1][0][alphas[i]]
                        }
                    }
                });
            }
        } else {
            //multiple predecessor
            let mut temp_alphas = HashMap::new();
            for pred in pred_hash.get(&i).unwrap() {
                let mut common_paths = paths_nodes[i].clone();
                common_paths.and(&paths_nodes[*pred]);

                //update alpha if needed
                if alphas[i] == paths_number + 1 {
                    if common_paths[alphas[*pred]] {
                        alphas[i] = alphas[*pred];
                    } else {
                        alphas[i] = common_paths.iter().position(|is_in| is_in).unwrap();
                    }
                }

                if common_paths[alphas[*pred]] {
                    let deltas = common_paths
                        .iter()
                        .enumerate()
                        .filter_map(|(path_id, is_in)| match is_in && path_id != alphas[*pred] {
                            true => Some(path_id),
                            false => None,
                        })
                        .collect::<Vec<usize>>();
                    temp_alphas.insert(alphas[*pred], deltas);

                    common_paths.iter().enumerate().for_each(|(path, is_in)| {
                        if is_in {
                            if path == alphas[*pred] {
                                dpm[i][0][alphas[*pred]] = dpm[*pred][0][alphas[*pred]]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                            } else {
                                dpm[i][0][path] = dpm[*pred][0][path];
                            }
                        }
                    });
                } else {
                    // new alpha needed
                    let temp_alpha = if common_paths[alphas[i]] {
                        alphas[i]
                    } else {
                        common_paths.iter().position(|is_in| is_in).unwrap()
                    };
                    let deltas = common_paths
                        .iter()
                        .enumerate()
                        .filter_map(|(path_id, is_in)| match is_in && path_id != alphas[*pred] {
                            true => Some(path_id),
                            false => None,
                        })
                        .collect::<Vec<usize>>();
                    temp_alphas.insert(alphas[*pred], deltas);

                    common_paths.iter().enumerate().for_each(|(path, is_in)| {
                        if is_in {
                            if path == temp_alpha {
                                dpm[i][0][temp_alpha] = dpm[*pred][0][alphas[*pred]]
                                    + dpm[*pred][0][temp_alpha]
                                    + score_matrix.get(&(lnz[i], '-')).unwrap();
                            } else {
                                dpm[i][0][path] = dpm[*pred][0][path] - dpm[*pred][0][temp_alpha];
                            }
                        }
                    });
                }
            }
            // remove temp alpha if needed

            if temp_alphas.keys().len() > 1 {
                for (a, delta) in temp_alphas.iter() {
                    if *a != alphas[i] {
                        dpm[i][0][*a] -= dpm[i][0][alphas[i]];
                        for path in delta.iter() {
                            if path != a {
                                dpm[i][0][*path] += dpm[i][0][*a];
                            }
                        }
                    }
                }
            }
        }
    }

    for j in 1..sequence.len() {
        dpm[0][j][alphas[0]] =
            dpm[0][j - 1][alphas[0]] + score_matrix.get(&(sequence[j], '-')).unwrap();
        for k in 0..paths_number {
            if k != alphas[0] {
                dpm[0][j][k] = dpm[0][j - 1][k];
            }
        }
    }
    //TODO: restart here
    for i in 1..lnz.len() - 1 {
        for j in 1..sequence.len() {}
    }
    /*
    dpm.iter().enumerate().for_each(|(i, line)| {
        print!("{}\tALPHA: {} {:?}\tDELTA: ",i, alphas[i],line[0][alphas[i]]);
        for p in 0..paths_number {
            if p != alphas[i] && paths_nodes[i][p] {
                print!("{} {:?}\t",p, line[0][alphas[i]]+line[0][p])
            }
        }
        println!()
    });

    println!("SCORE");
    println!("{}", dpm[lnz.len() - 2][0][alphas[lnz.len() - 2]]);
    for path in 0..paths_number {
        if path != alphas[lnz.len() - 2] {
            println!(
                "{}",
                dpm[lnz.len() - 2][0][alphas[lnz.len() - 2]] + dpm[lnz.len() - 2][0][path]
            );
        }
    }

    for i in 0..lnz.len() - 1 {
        for path in 0..paths_number {
            if paths_nodes[i][path] {
                if path == alphas[i] {
                    print!("{}\t", dpm[i][0][path])
                } else {
                    print!("{}\t", dpm[i][0][path] + dpm[i][0][alphas[i]])
                }
            } else {
                print!("x\t")
            }
        }
        println!()
    }
    */
}
