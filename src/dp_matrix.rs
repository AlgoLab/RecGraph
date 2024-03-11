use bit_vec::BitVec;

#[derive(Debug, Clone)]
pub struct DeltaScore {
    pub coord: Option<Box<[usize; 2]>>,
    pub values: Option<Box<[i32]>>,
}

impl Default for DeltaScore {
    fn default() -> Self {
        Self::new()
    }
}

impl DeltaScore {
    pub fn new() -> DeltaScore {
        DeltaScore {
            coord: None,
            values: None,
        }
    }
    pub fn new_coord(coord: (usize, usize)) -> DeltaScore {
        DeltaScore {
            coord: Some(Box::new([coord.0, coord.1])),
            values: None,
        }
    }

    pub fn new_values(values: Vec<i32>) -> DeltaScore {
        DeltaScore {
            coord: None,
            values: Some(values.into_boxed_slice()),
        }
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

    pub fn get_cell(&self, i: usize, j: usize) -> DeltaScore {
        self.deltas[i * self.cols + j].clone()
    }

    pub fn set_cell(&mut self, i: usize, j: usize, cell: DeltaScore) {
        self.deltas[i * self.cols + j] = cell;
    }

    pub fn get_val(&self, i: usize, j: usize, path: usize) -> i32 {
        let pred_cell = self.get_cell_ref(i, j);
        if pred_cell.values.as_ref().is_some() {
            pred_cell.values.as_ref().unwrap()[path]
        } else {
            let coord = pred_cell.coord.as_ref().unwrap();
            self.get_cell_ref(coord[0], coord[1])
                .values
                .as_ref()
                .unwrap()[path]
        }
    }

    pub fn get_cell_ref_mut(&mut self, i: usize, j: usize) -> &mut DeltaScore {
        &mut self.deltas[i * self.cols + j]
    }
    pub fn get_cell_ref(&self, i: usize, j: usize) -> &DeltaScore {
        &self.deltas[i * self.cols + j]
    }
    pub fn set_val(&mut self, i: usize, j: usize, path: usize, val: i32) {
        if self.get_cell_ref_mut(i, j).values.is_some() {
            self.get_cell_ref_mut(i, j).values.as_mut().unwrap()[path] = val
        } else {
            panic!("No value error")
        }
    }

    pub fn set_vals(&mut self, i: usize, j: usize, vals: Vec<i32>) {
        let cell = DeltaScore {
            coord: None,
            values: Some(vals.into_boxed_slice()),
        };
        self.set_cell(i, j, cell);
    }

    pub fn set_coord(&mut self, i: usize, j: usize, coord: (usize, usize)) {
        let cell = DeltaScore::new_coord(coord);
        self.set_cell(i, j, cell)
    }

    pub fn get_pred_coord(&self, i_pred: usize, j_pred: usize) -> (usize, usize) {
        let pred_cell = self.get_cell_ref(i_pred, j_pred);
        if pred_cell.coord.is_some() {
            let corrd = pred_cell.coord.as_ref().unwrap();
            (corrd[0], corrd[1])
        } else if pred_cell.values.is_some() {
            (i_pred, j_pred)
        } else {
            panic!("empty cell in {i_pred}, {j_pred}")
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
