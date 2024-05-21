use bstr::BString;
use priority_queue::DoublePriorityQueue;
use rayon::prelude::*;

pub fn get_minimizers(sequence: &BString, chunk_size: usize, m_len: usize) -> Vec<(usize, &[u8])> {
    let seeds: Vec<_> = sequence.chunks_exact(chunk_size).collect::<Vec<_>>();
    let minimizers: Vec<_> = seeds
        .par_iter()
        .map(|seed| {
            let mut hash_heap = DoublePriorityQueue::new();
            for i in 0..chunk_size - m_len + 1 {
                hash_heap.push(i, &seed[i..i + m_len]);
            }
            let (idx, minimizer) = hash_heap.peek_min().unwrap();
            (*idx, *minimizer)
        })
        .collect();
    minimizers
}
