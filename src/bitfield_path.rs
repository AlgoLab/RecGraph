use bitvec::prelude::*;

fn dir_u16_from_char(c: char) -> u16 {
    match c {
        'O' => 0,
        'D' => 1,
        'd' => 2,
        'L' => 3,
        'U' => 4,
        'X' => 5,
        'Y' => 6,
        'M' => 7,
        _ => panic! {"impossible direction char"},
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
        5 => 'X',
        6 => 'Y',
        7 => 'M',
        _ => panic! {"impossible direction bitslice"},
    }
}
pub fn pred_from_bitvec(bv: &BitVec<u16, Msb0>) -> usize {
    let pred: u16 = bv[..16].load_be();
    pred as usize
}

pub fn dir_from_bitvec(bv: &BitVec<u16, Msb0>) -> char {
    char_from_bitslice(&bv[16..])
}
pub fn set_path_cell(pred: usize, dir: char) -> BitVec<u16, Msb0> {
    let mut bv = bitvec![u16, Msb0; 0; 32];
    bv[..16].store::<u16>(pred as u16);
    bv[16..].store::<u16>(dir_u16_from_char(dir));
    bv
}
