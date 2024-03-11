/// GAFStruct represents a gaf alignment, with each field ordered as normal gaf field
#[derive(Debug, Clone)]
pub struct GAFStruct {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: char,
    pub path: Vec<usize>,
    pub path_length: usize,
    pub path_start: usize,
    pub path_end: usize,
    pub residue_matches_number: usize,
    pub alignment_block_length: String,
    pub mapping_quality: String,
    pub comments: String,
}
impl Default for GAFStruct {
    fn default() -> Self {
        Self::new()
    }
}

impl GAFStruct {
    pub fn new() -> GAFStruct {
        GAFStruct {
            query_name: String::from(""),
            query_length: 0,
            query_start: 0,
            query_end: 0,
            strand: ' ',
            path: vec![0usize],
            path_length: 0,
            path_start: 0,
            path_end: 0,
            residue_matches_number: 0,
            alignment_block_length: String::from(""),
            mapping_quality: String::from(""),
            comments: String::from(""),
        }
    }
    pub fn build_gaf_struct(
        query_name: String,
        query_length: usize,
        query_start: usize,
        query_end: usize,
        strand: char,
        path: Vec<usize>,
        path_length: usize,
        path_start: usize,
        path_end: usize,
        residue_matches_number: usize,
        alignment_block_length: String,
        mapping_quality: String,
        comments: String,
    ) -> GAFStruct {
        GAFStruct {
            query_name,
            query_length,
            query_start,
            query_end,
            strand,
            path,
            path_length,
            path_start,
            path_end,
            residue_matches_number,
            alignment_block_length,
            mapping_quality,
            comments,
        }
    }
    pub fn to_string(self) -> String {
        let path_matching: String = self
            .path
            .iter()
            .map(|id| id.to_string())
            .collect::<Vec<String>>()
            .join(">");
        let gaf_struct_to_string = format!(
            "{}\t{}\t{}\t{}\t{}\t>{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.query_name,
            self.query_length,
            self.query_start,
            self.query_end,
            self.strand,
            path_matching,
            self.path_length,
            self.path_start,
            self.path_end,
            self.residue_matches_number,
            self.alignment_block_length,
            self.mapping_quality,
            self.comments
        );
        gaf_struct_to_string
    }
}
