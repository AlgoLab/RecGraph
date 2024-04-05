use bit_vec::BitVec;

#[derive(Debug, Clone)]
pub enum DeltaScore {
    Coord(usize, usize),
    Values(Vec<i32>),
}

impl DeltaScore {
    pub fn get_val(&self, path: usize) -> i32 {
        match self {
            DeltaScore::Values(values) => values[path],
            DeltaScore::Coord(_, _) => panic!("No value error"),
        }
    }

    pub fn set_val(&mut self, path: usize, val: i32) {
        if let DeltaScore::Values(values) = self {
            values[path] = val;
        } else {
            panic!("No value error")
        }
    }

    pub fn new() -> DeltaScore {
        DeltaScore::Values(Vec::new())
    }
}
#[derive(Debug, Clone)]
pub struct DpMatrix {
    pub rows: usize,
    pub cols: usize,
    pub dpm: Vec<i32>,
}

impl DpMatrix {
    pub fn new(rows: usize, cols: usize) -> DpMatrix {
        DpMatrix {
            rows,
            cols,
            dpm: vec![0; rows * cols],
        }
    }

    pub fn empty_new() -> DpMatrix {
        DpMatrix {
            rows: 0,
            cols: 0,
            dpm: Vec::new(),
        }
    }

    pub fn get(&self, i: usize, j: usize) -> i32 {
        self.dpm[i * self.cols + j]
    }

    pub fn set(&mut self, i: usize, j: usize, value: i32) {
        self.dpm[i * self.cols + j] = value;
    }
}

#[derive(Debug, Clone)]
pub struct DpDeltas {
    pub rows: usize,
    pub cols: usize,
    pub deltas: Vec<DeltaScore>,
}

impl DpDeltas {
    pub fn new(rows: usize, cols: usize) -> DpDeltas {
        DpDeltas {
            rows,
            cols,
            deltas: vec![DeltaScore::new(); rows * cols],
        }
    }

    pub fn empty_new() -> DpDeltas {
        DpDeltas {
            rows: 0,
            cols: 0,
            deltas: Vec::new(),
        }
    }

    pub fn get_cell(&self, i: usize, j: usize) -> &DeltaScore {
        &self.deltas[i * self.cols + j]
    }

    pub fn set_cell(&mut self, i: usize, j: usize, cell: DeltaScore) {
        self.deltas[i * self.cols + j] = cell;
    }

    pub fn get_val(&self, i: usize, j: usize, path: usize) -> i32 {
        match self.get_cell(i, j) {
            DeltaScore::Values(values) => values[path],
            DeltaScore::Coord(x, y) => {
                let pred_cell = self.get_cell(*x, *y);
                if let DeltaScore::Values(vals) = pred_cell {
                    vals[path]
                } else {
                    panic!("No value error")
                }
            }
        }
    }

    pub fn set_val(&mut self, i: usize, j: usize, path: usize, val: i32) {
        if let DeltaScore::Values(values) = &mut self.deltas[i * self.cols + j] {
            values[path] = val;
        } else {
            panic!("No value error")
        }
    }

    pub fn set_vals(&mut self, i: usize, j: usize, vals: Vec<i32>) {
        self.deltas[i * self.cols + j] = DeltaScore::Values(vals);
    }

    pub fn set_coord(&mut self, i: usize, j: usize, coord: (usize, usize)) {
        self.deltas[i * self.cols + j] = DeltaScore::Coord(coord.0, coord.1);
    }

    pub fn get_pred_coord(&self, i_pred: usize, j_pred: usize) -> (usize, usize) {
        match self.get_cell(i_pred, j_pred) {
            DeltaScore::Coord(x, y) => (*x, *y),
            DeltaScore::Values(_) => (i_pred, j_pred),
        }
    }
}

pub fn print_abs_score_matrix(
    matrix: &DpMatrix,
    deltas: &DpDeltas,
    alphas: &Vec<usize>,
    node_paths: &Vec<BitVec>,
) {
    let paths_num = node_paths[0].len();
    for i in 1..matrix.rows - 1 {
        for j in 0..matrix.cols {
            print!("[");
            for path in 0..paths_num {
                if node_paths[i][path] {
                    if path == alphas[i] {
                        print!("{} ", matrix.get(i, j));
                    } else {
                        print!("{} ", deltas.get_val(i, j, path) + matrix.get(i, j));
                    }
                } else {
                    print!("x ");
                }
            }
            print!("] ");
        }
        println!();
    }
    println!();
}
