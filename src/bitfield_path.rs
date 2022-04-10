use bitfield::bitfield;

use crate::bitfield_path;


bitfield! {
    struct BestChoice (MSB0 [u32]);
    u32;
    get_pred, set_pred: 16, 0;
    get_dir, set_dir: 21, 17;
}

fn dir(c: char) -> u32 {
    match c {
        'D' => {
            0
        }
        'd' => {
            1
        }
        'L' => {
            2
        }
        'U' => {
            3
        }
        'O' => {
            4
        }
        _ => {
            5
        }
    }
}

pub fn test1() {
    let mut best_choice = BestChoice([0, 0]);
    best_choice.set_pred(1000);
    best_choice.set_dir(dir('d'));

    assert_eq!(best_choice.get_pred(),1000);
    assert_eq!(best_choice.get_dir(),1);
}

pub fn example_path() {
    let mut path: Vec<Vec<BestChoice<[u32; 2]>>> = Vec::with_capacity(2);
    let mut path_line = Vec::with_capacity(3);
    path_line.insert(0, BestChoice([0, dir('O')]));
    path_line.insert(1, BestChoice([0, dir('L')]));
    path_line.insert(2, BestChoice([0, dir('L')]));
    path.insert(0,path_line);
    
    let mut path_line = Vec::with_capacity(3);
    path_line.insert(0, BestChoice([0, dir('U')]));
    path_line.insert(1, BestChoice([0, dir('D')]));
    path_line.insert(2, BestChoice([1, dir('L')]));
    path.insert(1, path_line);
}

fn build_path(path: Vec<Vec<BestChoice<[u32; 2]>>>) {
    let mut row = path.len() - 1;
    let mut col = path[row].len() - 1;
    let mut moves = Vec::new();
    while path[row][col].get_dir() != 4 {
        match path[row][col].get_dir() {
            0 | 1=> {
                moves.push('D');
                row -= 1;
                col -= 1;
            }
            2 => {
                moves.push('L');
                col -= 1;
            }
            3 => {
                moves.push('U');
                row -= 1;
            }
            _ => {
                panic!();
            }
        }
    }
    assert_eq!('L', moves[0]);
    assert_eq!('D', moves[1]);
}


