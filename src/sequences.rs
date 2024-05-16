use bstr::BString;
use needletail::{parse_fastx_file, Sequence};

pub fn get_sequences(file_path: String) -> (Vec<BString>, Vec<BString>) {
    let mut sequences = Vec::new();
    let mut ids = Vec::new();

    let current_absolute_position = std::env::current_dir().unwrap();
    let mut reader = parse_fastx_file(&file_path).expect(&format!(
        "Invalid sequence path {} from {}",
        file_path,
        current_absolute_position.display()
    ));
    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid sequence");
        let seqrec_norm = seqrec.normalize(false);

        let mut sequence = BString::from(seqrec_norm.sequence());
        sequence.insert(0, b'$');
        sequence.push(b'$');
        sequences.push(sequence);
        ids.push(BString::from(seqrec.id()));
    }
    (sequences, ids)
}
