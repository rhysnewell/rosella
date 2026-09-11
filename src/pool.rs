use std::sync::{Arc, OnceLock};

use anyhow::Result;
use frugal::api::META_PREDICTOR_STACK_SIZE;
use rayon::ThreadPool;

static POOL: OnceLock<Arc<ThreadPool>> = OnceLock::new();

/// One pool for the whole process. The gene caller needs an owned pool it can be handed, which
/// rayon's global pool cannot supply, so nothing builds a second one and everything installs.
pub fn init(threads: usize) -> Result<()> {
    let built = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .stack_size(META_PREDICTOR_STACK_SIZE)
        .build()?;
    let _ = POOL.set(Arc::new(built));
    Ok(())
}

pub fn get() -> Arc<ThreadPool> {
    POOL.get()
        .expect("the thread pool is built before any subcommand runs")
        .clone()
}

pub fn install<R: Send>(op: impl FnOnce() -> R + Send) -> R {
    get().install(op)
}

pub fn threads() -> usize {
    get().current_num_threads()
}
