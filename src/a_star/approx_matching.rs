use bio::pattern_matching::ukkonen;

use lt_fm_index::LtFmIndexBuilder;


pub fn test(seq: &Vec<u8>, pattern: &Vec<u8>) {

    let start = std::time::Instant::now();
    let mut uk = ukkonen::Ukkonen::with_capacity(10, ukkonen::unit_cost);
    let occ: Vec<(usize, usize)> = uk.find_all_end(pattern, seq, 1).collect();
    println!("Ukkonen {:?}", start.elapsed());

    let seq = seq.clone();
    let start = std::time::Instant::now();
     let builder = LtFmIndexBuilder::new()
                .text_type_is_inferred()
                .set_lookup_table_kmer_size_to_default()
                .set_suffix_array_sampling_ratio_to_default();    
    
    let fm_index = builder.build(seq).unwrap();
    let matches_pos = fm_index.locate(&pattern);
    println!("FM {:?}", start.elapsed());
}
