use bstr::BString;
use needletail::{parse_fastx_file, Sequence};

pub fn get_sequences(file_path: String) -> (Vec<BString>, Vec<BString>) {
    let mut sequences = Vec::new();
    let mut ids = Vec::new();

    let mut reader = parse_fastx_file(&file_path).expect("Invalid sequence path");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid sequence");
        let seqrec_norm = seqrec.normalize(false);

        let mut sequence = BString::from(seqrec_norm.sequence());
        sequence.insert(0, b'$');
        sequences.push(sequence);
        ids.push(BString::from(seqrec.id()));
    }
    (sequences, ids)
}
