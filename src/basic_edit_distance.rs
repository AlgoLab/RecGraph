use std::cmp;

use ndarray::Array2;

pub fn exec_ed(s1: &Vec<char>, s2: &Vec<char>) {
    let (dp_mat, _pred_mat) = pd_ed(&s1, &s2);
    
    
    println!("EDIT DISTANCE: {}", dp_mat[[s1.len()-1,s2.len()-1]]);
    println!("");
}

fn pd_ed(s1: &Vec<char>, s2: &Vec<char>) -> (Array2<i32>, Vec<Vec<char>>) {
    let mut a = Array2::<i32>::zeros((s1.len(), s2.len()));
    let mut path =vec![vec!['x';s2.len()];s1.len()];
    
    for (m, char_at_m) in s1.iter().enumerate() {
        for (n, char_at_n) in s2.iter().enumerate() {
            match (m, n) {
                (0, 0) => {
                    a[[0,0]] = 0 as i32;
                    path [0][0] = 'x';
                }
                (_, 0) => {
                    a[[m,0]] = m as i32;
                    path [m][0] = 'U';
                },
                (0, _) => {
                    a[[0,n]] = n as i32;
                    path [0][n] = 'L';
                },
                _ => {
                    let mut to_add = 0;
                    if char_at_m != char_at_n {
                    to_add = 1;
                    }
            
                    a[[m,n]] = cmp::min(a[[m-1,n-1]]+to_add, cmp::min(a[[m,n-1]]+1, a[[m-1,n]]+1));
                    match a[[m,n]] {
                        _ if a[[m,n]] == a[[m-1,n-1]] && char_at_m == char_at_n => path[m][n] = 'D',
                        _ if a[[m,n]] == a[[m-1,n-1]] + 1 => path[m][n] = 'd',
                        _ if a[[m,n]] == a[[m,n-1]] + 1 => path[m][n] = 'L',
                        _ if a[[m,n]] == a[[m-1,n]] + 1 => path[m][n] = 'U',
                        _ => panic!("impossible")
                    }
                    
                }
            }
        }
    }
    (a, path)
}
