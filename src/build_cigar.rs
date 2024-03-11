pub fn build_cigar(cigar: &Vec<char>) -> String {
    let mut output_string = String::new();

    let mut d_count = 0;
    let mut u_count = 0;
    let mut l_count = 0;
    let mut mm_count = 0;
    for ch in cigar.iter() {
        match ch {
            'D' => {
                if u_count != 0 {
                    output_string = format!("{}{}I", output_string, u_count);
                    u_count = 0
                }
                if l_count != 0 {
                    output_string = format!("{}{}D", output_string, l_count);
                    l_count = 0
                }
                if mm_count != 0 {
                    output_string = format!("{}{}X", output_string, mm_count);
                    mm_count = 0
                }
                d_count += 1;
            }

            'U' => {
                if d_count != 0 {
                    output_string = format!("{}{}M", output_string, d_count);
                    d_count = 0
                }
                if l_count != 0 {
                    output_string = format!("{}{}D", output_string, l_count);
                    l_count = 0
                }
                if mm_count != 0 {
                    output_string = format!("{}{}X", output_string, mm_count);
                    mm_count = 0
                }
                u_count += 1;
            }
            'd' => {
                if d_count != 0 {
                    output_string = format!("{}{}M", output_string, d_count);
                    d_count = 0
                }
                if l_count != 0 {
                    output_string = format!("{}{}D", output_string, l_count);
                    l_count = 0
                }
                if u_count != 0 {
                    output_string = format!("{}{}I", output_string, u_count);
                    u_count = 0
                }
                mm_count += 1;
            }
            _ => {
                if d_count != 0 {
                    output_string = format!("{}{}M", output_string, d_count);
                    d_count = 0;
                }
                if u_count != 0 {
                    output_string = format!("{}{}I", output_string, u_count);
                    u_count = 0
                }
                if mm_count != 0 {
                    output_string = format!("{}{}X", output_string, mm_count);
                    mm_count = 0
                }
                l_count += 1;
            }
        }
    }
    if d_count != 0 {
        output_string = format!("{}{}M", output_string, d_count);
    }
    if u_count != 0 {
        output_string = format!("{}{}I", output_string, u_count);
    }
    if l_count != 0 {
        output_string = format!("{}{}D", output_string, l_count);
    }
    if mm_count != 0 {
        output_string = format!("{}{}X", output_string, mm_count);
    }
    output_string
}

