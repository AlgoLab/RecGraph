use bstr::BString;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use handlegraph::hashgraph::HashGraph;
use lt_fm_index::LtFmIndex;
use lt_fm_index::blocks::Block3;
use needletail::Sequence;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::io::{self, Error, Read};
use std::time::{Duration, Instant};

use crate::a_star::{a_star_output, a_star_visit, chain_heur, fast_heuristic, seed_heuristic};
use crate::args_parser::{ClArgs, EstimateFunction, ScoringParams};
use crate::new_path_graph::path_graph::{PathGraph, remove_duplicate_paths};
use crate::sequences;

use super::a_star_output::build_gaf;

pub fn a_star_demo_chain() -> Result<(), Error> {
    let args = ClArgs::parse();
    configure_rayon_threads(args.threads)?;
    let file_path = args.graph_path;
    if file_path.is_empty() {
        return Err(Error::new(
            std::io::ErrorKind::InvalidInput,
            "Graph path is empty",
        ));
    }
    let mut stdin_used = false;
    let (gfa_data, sequence_data) = if file_path == "-" && args.sequence_path == "-" {
        let (gfa_data, sequence_data) = get_input_from_stdin_new()?;
        stdin_used = true;
        (gfa_data, sequence_data)
    } else {
        (GFA::new(), (Vec::new(), Vec::new()))
    };

    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = if file_path == "-" && !stdin_used {
        let stdin_buf = StdinLineBuffer::new()?;
        parser
            .parse_lines(stdin_buf.iter())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
    } else if !stdin_used {
        parser.parse_file(file_path).unwrap()
    } else {
        gfa_data
    };
    eprintln!(
        "Parsed GFA with {} segments and {} paths",
        gfa.segments.len(),
        gfa.paths.len()
    );

    let mut graph: HashGraph = HashGraph::from_gfa(&gfa);
    let original_path_ids = remove_duplicate_paths(&mut graph);
    let path_graph = PathGraph::from_hash_graph(&graph);
    let (sequences, names) = if !stdin_used {
        sequences::get_sequences(args.sequence_path.clone(), args.sequence_path == "-")
    } else {
        sequence_data
    };
    let chunk_size = args.seed_len;
    let indexes = path_graph.get_indexes();
    let mut outs = Vec::new();
    //let mut explored_pos_vec = Vec::new();
    let mut heur_tot_time = Duration::new(0, 0);
    let mut explore_tot_time = Duration::new(0, 0);
    sequences.iter().zip(names).for_each(|(seq, name)| {
        let istant = Instant::now();
        let mut heuristic = match args.est_function {
            EstimateFunction::Chaining => chain_heur::build_heuristic(
                &indexes,
                seq,
                chunk_size as usize,
                args.base_rec_cost as usize,
                &path_graph,
                args.max_rec > 0,
            ),
            EstimateFunction::Seeding => seed_heuristic::build_heuristic(
                &indexes,
                seq,
                chunk_size as usize,
                args.base_rec_cost as usize,
                &path_graph,
                args.max_rec > 0,
            ),
            EstimateFunction::Fast => fast_heuristic::build_heuristic(
                &indexes,
                seq,
                chunk_size as usize,
                args.base_rec_cost as usize,
                &path_graph,
                args.max_rec > 0,
            ),
        };

        let matches_in_path: Vec<_> = heuristic
            .1
            .par_iter()
            .map(|m| a_star_output::get_matches_end_in_path(m))
            .collect();
        let heu_time = istant.elapsed();
        heur_tot_time += heu_time;

        let explore_start = Instant::now();
        let (end_pos, mut alignment_graph /* , explored_pos*/) = a_star_visit::exec(
            seq,
            &mut heuristic,
            &path_graph,
            args.alignment_mode,
            args.base_rec_cost as u16,
            args.max_rec as u32,
            &args.est_function,
            &args.scoring_params,
        );
        let explore_time = explore_start.elapsed();
        explore_tot_time += explore_time;

        //explored_pos_vec.push(explored_pos);
        let graph_size = (path_graph.get_graph_size() * seq.len()) as f32;
        let explored_cells = (alignment_graph.len() as f32 / graph_size) * 100.0;

        outs.push(build_gaf(
            &mut alignment_graph,
            &end_pos,
            &path_graph,
            seq,
            args.alignment_mode,
            &matches_in_path,
            chunk_size as usize,
            &indexes,
            &name,
            explore_time + heu_time,
            &original_path_ids,
            explored_cells,
        ));
    });
    //a_star_output::save_coords(&args.out_file, &explored_pos_vec);
    outs.iter().for_each(|out| println!("{}", out.to_string()));
    let mem = peak_mem_usage().unwrap();

    eprintln!("Peak memory (Byte)\t{}", mem);
    eprintln!("Heuristic tot time\t{:?}", heur_tot_time);
    eprintln!("Explore tot time\t{:?}", explore_tot_time);
    /*
    if args.alignment_mode {
        crate::a_star::check_ed::semiglobal_test(&path_graph, &sequences)
    } else {
        crate::a_star::check_ed::test(&path_graph, &sequences)
    }
    */
    Ok(())
}

fn configure_rayon_threads(threads: Option<usize>) -> Result<(), Error> {
    if let Some(threads) = threads {
        if threads == 0 {
            return Err(Error::new(
                io::ErrorKind::InvalidInput,
                "--threads must be greater than 0",
            ));
        }

        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .map_err(|e| {
                Error::new(
                    io::ErrorKind::Other,
                    format!("Failed to configure Rayon thread pool: {e}"),
                )
            })?;
    }
    Ok(())
}

#[cfg(target_os = "linux")]
fn peak_mem_usage() -> Result<usize, &'static str> {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        let retval = libc::getrusage(libc::RUSAGE_SELF, &mut rusage as *mut _);
        match retval {
            0 => Ok(rusage.ru_maxrss as usize * 1024),
            _ => Err("Error getting memory usage"),
        }
    }
}

pub fn alignment_bench(
    sequences: &Vec<BString>,
    path_graph: &PathGraph,
    chunk_size: u32,
    rec_cost: usize,
    is_local: bool,
    max_rec: u32,
    indexes: &Vec<(LtFmIndex<u32, Block3<u128>>, Vec<u32>)>,
    est_function: EstimateFunction,
) {
    sequences.iter().for_each(|seq| {
        let mut heuristic = match est_function {
            EstimateFunction::Chaining => chain_heur::build_heuristic(
                indexes,
                seq,
                chunk_size as usize,
                rec_cost,
                path_graph,
                max_rec > 0,
            ),
            EstimateFunction::Seeding => seed_heuristic::build_heuristic(
                indexes,
                seq,
                chunk_size as usize,
                rec_cost,
                path_graph,
                max_rec > 0,
            ),
            EstimateFunction::Fast => fast_heuristic::build_heuristic(
                indexes,
                seq,
                chunk_size as usize,
                rec_cost,
                path_graph,
                max_rec > 0,
            ),
        };
        a_star_visit::exec(
            seq,
            &mut heuristic,
            &path_graph,
            is_local,
            rec_cost as u16,
            max_rec,
            &est_function,
            &ScoringParams::base(),
        );
    });
}

pub struct StdinLineBuffer {
    buffer: Vec<u8>,
}

impl StdinLineBuffer {
    pub fn new() -> io::Result<Self> {
        let mut buffer = Vec::new();
        io::stdin().read_to_end(&mut buffer)?;
        Ok(Self { buffer })
    }

    pub fn iter(&self) -> impl Iterator<Item = &[u8]> {
        self.buffer.split(|&b| b == b'\n').map(|line| {
            if line.ends_with(b"\r") {
                &line[..line.len() - 1]
            } else {
                line
            }
        })
    }
}

fn get_input_from_stdin() -> io::Result<(GFA<usize, ()>, (Vec<BString>, Vec<BString>))> {
    // Keep reading until line starts with '>', which indicates start of sequences and end of GFA

    let mut gfa_data = String::new();
    let mut sequence_data = String::new();
    for line in io::stdin().lines() {
        let line = line?;
        if line.starts_with('>') {
            sequence_data.push_str(&line);
            sequence_data.push('\n');
            break;
        }
        gfa_data.push_str(&line);
        gfa_data.push('\n');
    }
    for line in io::stdin().lines() {
        let line = line?;
        sequence_data.push_str(&line);
        sequence_data.push('\n');
    }
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser
        .parse_lines(gfa_data.lines())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut seq_reader = needletail::parse_fastx_reader(sequence_data.as_bytes())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut sequences = Vec::new();
    let mut ids = Vec::new();

    while let Some(record) = seq_reader.next() {
        let seqrec = record.expect("Invalid sequence");
        let seqrec_norm = seqrec.normalize(true);

        let mut sequence = BString::from(seqrec_norm.sequence());
        sequence.insert(0, b'$');
        sequence.push(b'$');
        sequences.push(sequence);
        ids.push(BString::from(seqrec.id()));
    }
    Ok((gfa, (sequences, ids)))
}

fn get_input_from_stdin_new() -> io::Result<(GFA<usize, ()>, (Vec<BString>, Vec<BString>))> {
    // Keep reading until line starts with '>', which indicates start of sequences and end of GFA

    let mut gfa_data = String::new();
    let mut sequence_data = String::new();
    let mut lines = io::stdin().lines();
    if let Some(header) = lines.next() {
        let header = header?;
        if header != "GRAPH:" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Wrong input format: expected 'GRAPH:' line",
            ));
        } else {
            while let Some(line) = lines.next() {
                let line = line?;
                if line == "READ:" {
                    break;
                }
                gfa_data.push_str(&line);
                gfa_data.push('\n');
            }
            while let Some(line) = lines.next() {
                let line = line?;
                sequence_data.push_str(&line);
                sequence_data.push('\n');
            }
        }
    } else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No input provided",
        ));
    }
    
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser
        .parse_lines(gfa_data.lines())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut seq_reader = needletail::parse_fastx_reader(sequence_data.as_bytes())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let mut sequences = Vec::new();
    let mut ids = Vec::new();

    while let Some(record) = seq_reader.next() {
        let seqrec = record.expect("Invalid sequence");
        let seqrec_norm = seqrec.normalize(true);

        let mut sequence = BString::from(seqrec_norm.sequence());
        sequence.insert(0, b'$');
        sequence.push(b'$');
        sequences.push(sequence);
        ids.push(BString::from(seqrec.id()));
    }
    Ok((gfa, (sequences, ids)))
}
