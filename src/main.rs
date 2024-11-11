use recalign::a_star::a_star_demo;
use std::time::SystemTime;

#[cfg(target_os = "linux")]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

fn main() {
    let now = SystemTime::now();

    a_star_demo::a_star_demo_chain();
    
    match now.elapsed() {
        Ok(elapsed) => {
            // it prints '2'
            eprintln!("Done in {:?}.", elapsed);
        }
        Err(e) => {
            // an error occurred!
            eprintln!("Error: {e:?}");
        }
    }
}
