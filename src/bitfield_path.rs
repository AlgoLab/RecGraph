use bitvec::{prelude::*, view::BitView};

pub fn test1() {
    let mut bv = bitvec![u16, Msb0; 0; 32];
    bv[..16].store::<u16>(1821); //first 16 bits pred
    bv[16..].store::<u16>(dir_u16_from_char('L')); //last 16 bits (will be last 5) direction
    let (pred, dir) = bv.split_at(16);
    let pred: u16 = pred.load_be();
    assert_eq!(pred, 1821);
    assert_eq!('L', char_from_bitslice(dir));

}

pub fn example_path() {
    let bv = bitvec![u16, Msb0; 0; 32];
    let mut path = vec![vec![bv; 5]; 5];
    for i in 1..path.len() {
        path[i][0][..16].store::<u16>((i-1) as u16);
        path[i][0][16..].store::<u16>(dir_u16_from_char('U'));
    }
    for j in 1..path[0].len() {
        path[0][j][16..].store::<u16>(dir_u16_from_char('L'));
    }

    for i in 1..path.len(){
        for j in 1..path[i].len() {
            path[i][j][..16].store::<u16>((i-1) as u16);
            if (i+j)%3 == 0 {
                path[i][j][16..].store::<u16>(dir_u16_from_char('D'));
            } else {
                path[i][j][16..].store::<u16>(dir_u16_from_char('d'));
            }
        }
    }
    build_path(&path);
}
fn build_path(path: &Vec<Vec<BitVec<u16, Msb0>>>) {
    let mut row = path.len() - 1;
    let mut col = path[row].len() - 1;
    let mut moves = Vec::new();
    while char_from_bitslice(&path[row][col][16..]) != 'O' {
        let dir = char_from_bitslice(&path[row][col][16..]);
        moves.push(dir);
        match dir {
            'D' | 'd' => {
                col -= 1;
                let pred: u16 = path[row][col][..16].load_be();
                row = pred as usize;
            }
            'U' => {
                let pred: u16 = path[row][col][..16].load_be();
                row = pred as usize;
            }
            'L' => {
                col -= 1;
            }
            _ => panic!("impossible value while building path")
        }
    }
    let expected = ['d', 'D', 'd', 'd'];
    for i in 0..moves.len() {
        assert_eq!(moves[i], expected[i])
    }
}
fn dir_u16_from_char(c: char) -> u16 {
    match c {
        'O' => 0,
        'D' => 1,
        'd' => 2,
        'L' => 3,
        'U' => 4,
        _  => panic!{"impossible direction char"}
    }
}

fn char_from_bitslice(bs: &BitSlice<u16, Msb0>) -> char {
    let dir: u16 = bs.load_be();
    match dir {
        0 => 'O',
        1 => 'D',
        2 => 'd',
        3 => 'L',
        4 => 'U',
        _ => panic!{"impossible direction bitslice"}
    }
}