use std::cmp;

pub fn ab_ed(s1: &Vec<char>, s2: &Vec<char>, ampl: i32, a: &mut Vec<Vec<i32>>, path: &mut Vec<Vec<char>>) {
    for i in 0..s1.len() as i32 {
        for j in i-ampl..=i+ampl{
            let i_us = i as usize;
            let j_us = j as usize;
            match (i,j){
                _ if j<0 || j>=s2.len() as i32 => {},
                (0,0) => {
                    path[0][0] = 'O'
                },
                (_, 0) => {
                    a[i_us][0] = i;
                    path [i_us][0] = 'U';
                },
                (0, _) => {
                    a[0][j_us] = j;
                    path [0][j_us] = 'L';
                },
                _ if i-j==ampl => {
                    let mut to_add = 0;
                    if s1[i_us] != s2[j_us] {
                    to_add = 1;
                    }
            
                    a[i_us][j_us] = cmp::min(a[i_us-1][j_us-1]+to_add, 
                                             a[i_us-1][j_us]+1);
                    match a[i_us][j_us] {
                        _ if a[i_us][j_us] == a[i_us-1][j_us-1] && s1[i_us] == s2[j_us] => path[i_us][j_us] = 'D',
                        _ if a[i_us][j_us] == a[i_us-1][j_us-1] + 1 => path[i_us][j_us] = 'd',
                        _ if a[i_us][j_us] == a[i_us-1][j_us] + 1 => path[i_us][j_us] = 'U',
                        _ => panic!("impossible")
                    }
                },
                _ if i-j==-ampl => {
                    let mut to_add = 0;
                    if s1[i_us] != s2[j_us] {
                    to_add = 1;
                    }
            
                    a[i_us][j_us] = cmp::min(a[i_us-1][j_us-1]+to_add, 
                                             a[i_us][j_us-1]+1);
                    match a[i_us][j_us] {
                        _ if a[i_us][j_us] == a[i_us-1][j_us-1] && s1[i_us] == s2[j_us] => path[i_us][j_us] = 'D',
                        _ if a[i_us][j_us] == a[i_us-1][j_us-1] + 1 => path[i_us][j_us] = 'd',
                        _ if a[i_us][j_us] == a[i_us][j_us-1] + 1 => path[i_us][j_us] = 'L',
                        _ => panic!("impossible")
                    }
                },
                _ => {
                    let mut to_add = 0;
                    if s1[i_us] != s2[j_us] {
                    to_add = 1;
                    }
            
                    a[i_us][j_us] = cmp::min(a[i_us-1][j_us-1]+to_add, 
                                             cmp::min(a[i_us-1][j_us]+1, a[i_us][j_us-1]+1));
                    match a[i_us][j_us] {
                        _ if a[i_us][j_us] == a[i_us-1][j_us-1] && s1[i_us] == s2[j_us] => path[i_us][j_us] = 'D',
                        _ if a[i_us][j_us] == a[i_us-1][j_us-1] + 1 => path[i_us][j_us] = 'd',
                        _ if a[i_us][j_us] == a[i_us-1][j_us] + 1 => path[i_us][j_us] = 'U',
                        _ if a[i_us][j_us] == a[i_us][j_us-1] + 1 => path[i_us][j_us] = 'L',

                        _ => panic!("impossible")
                    }
                },
            }
        } 
    }
    match ampl_is_enough_iterative(path) {
        true => {
            println!("EDIT DISTANCE: {}", a[s1.len()-1][s2.len()-1]);},
        false => {
            ab_ed(s1, s2, ampl*2, a, path)}
    }    
}

fn ampl_is_enough_iterative(path: &Vec<Vec<char>>) -> bool{
    let mut row = path.len()-1;
    let mut col = path[0].len()-1;
    while path[row][col] != 'O' {
        if row == 0 || col == 0 {return true}
        else if path[row-1][col-1] == 'x' ||
           path[row][col-1] == 'x' ||
           path[row-1][col] == 'x' {
                   return false
               } 
        else {
            match path[row][col] {
                'D'|'d' => {
                    row -= 1;
                    col -= 1;
                },
                'U' => {
                    row -= 1;
                },
                'L' => {
                    col -= 1;
                },
                _ => panic!("ampl_is_enough panic")
            }
               }
    }
    return true
}