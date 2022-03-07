use std::cmp;

pub fn ab_ed_km(s1: &[char], s2: &[char], ampl: i32) {
    let s1_len = s1.len();
    let s2_len = s2.len();
    let mut a = vec![vec![(s1_len + s2_len) as i32; ampl as usize]; s1_len];
    let mut path = vec![vec!['x'; ampl as usize]; s1_len];

    for i in 0..s1_len {
        for j in 0..ampl as usize {
            let i_32 = i as i32;
            let j_32 = j as i32;
            match (i, j) {
                _ if i_32 + j_32 - ampl / 2 >= s2_len as i32 || i_32 + j_32 < ampl / 2 => {} // celle angolo sup sx e inf dx fuori dalla banda
                _ if i == 0 => {
                    if i_32 + j_32 != ampl / 2 {
                        // prima riga tabella
                        a[0][j] = a[0][j - 1] + 1;
                        path[0][j] = 'L';
                    } else {
                        // origine
                        a[i][j] = 0;
                        path[i][j] = 'O'
                    }
                }
                (_, 0) => {
                    // bordo inf banda
                    let mut to_add = 0;
                    if s1[i] != s2[i + j - (ampl / 2) as usize] {
                        to_add = 1
                    }
                    a[i][0] = cmp::min(a[i - 1][0] + to_add, a[i - 1][j + 1] + 1);
                    if a[i][0] == a[i - 1][0] && to_add == 0 {
                        path[i][0] = 'D'
                    } else if a[i][0] == a[i - 1][0] + 1 {
                        path[i][0] = 'd'
                    } else {
                        path[i][0] = 'U';
                    }
                }
                _ if j == ampl as usize - 1 => {
                    //bordo sup banda
                    let mut to_add = 0;
                    if s1[i] != s2[i + j - (ampl / 2) as usize] {
                        to_add = 1
                    }
                    a[i][j] = cmp::min(a[i - 1][j] + to_add, a[i][j - 1] + 1);
                    if a[i][j] == a[i - 1][j] && to_add == 0 {
                        path[i][j] = 'D'
                    } else if a[i][j] == a[i - 1][j] + 1 {
                        path[i][j] = 'd'
                    } else {
                        path[i][j] = 'L';
                    }
                }
                _ => {
                    // celle interne alla banda
                    let mut to_add = 0;
                    if s1[i] != s2[i + j - (ampl / 2) as usize] {
                        to_add = 1
                    }
                    a[i][j] = cmp::min(
                        a[i - 1][j] + to_add,
                        cmp::min(a[i][j - 1] + 1, a[i - 1][j + 1] + 1),
                    );
                    if a[i][j] == a[i - 1][j] && to_add == 0 {
                        path[i][j] = 'D'
                    } else if a[i][j] == a[i - 1][j] + 1 {
                        path[i][j] = 'd'
                    } else if a[i][j] == a[i][j - 1] + 1 {
                        path[i][j] = 'L';
                    } else {
                        path[i][j] = 'U';
                    }
                }
            }
        }
    }
    match ampl_is_enough_iterative(&path, s2_len - 1 + (ampl / 2) as usize - (s1_len - 1)) {
        true => {
            println!(
                "EDIT DISTANCE: {}",
                a[s1_len - 1][s2_len - 1 + (ampl / 2) as usize - (s1_len - 1)]
            );
        }
        false => ab_ed_km(s1, s2, ampl * 2 + 1),
    }
}

fn ampl_is_enough_iterative(path: &[Vec<char>], start_col: usize) -> bool {
    let mut row = path.len() - 1;
    let mut col = start_col;
    let col_number = path[0].len();

    while path[row][col] != 'O' {
        if col == 0 || col == col_number - 1 {
            if path[row][col] == 'D' {
                row -= 1; // finchè ho match posso continuare anche se sul margine
            } else {
                return false;
            }
        } else {
            match path[row][col] {
                'D' | 'd' => {
                    row -= 1;
                }
                'U' => {
                    row -= 1;
                    col += 1;
                }
                'L' => {
                    col -= 1;
                }
                _ => panic!("ampl_is_enough panic"),
            }
        }
    }
    true
}
