use std::ffi::{CString, OsStr, OsString, c_char, c_int};
use std::io::{BufReader, PipeReader};
use std::os::fd::IntoRawFd;
use std::os::unix::ffi::OsStrExt;
use std::sync::Mutex;

use anyhow::{Result, bail};
// Only the C library calls zlib, and rustc leaves out a crate no Rust code names.
use libz_sys as _;

unsafe extern "C" {
    fn rosella_minimap2(argc: c_int, argv: *mut *mut c_char, out: c_int) -> c_int;
}

// The library keeps its output stream and log level in globals, so one run at a time.
static RUNNING: Mutex<()> = Mutex::new(());

// Its log would print over the progress bars, so only its errors are kept.
pub fn minimap2<T: Send>(
    args: &[OsString],
    read: impl FnOnce(&mut BufReader<PipeReader>) -> Result<T> + Send,
) -> Result<T> {
    let owned = [OsStr::new("minimap2"), OsStr::new("-v"), OsStr::new("1")]
        .into_iter()
        .chain(args.iter().map(OsString::as_os_str))
        .map(|arg| CString::new(arg.as_bytes()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut argv: Vec<*mut c_char> = owned.iter().map(|arg| arg.as_ptr().cast_mut()).collect();
    let argc = argv.len() as c_int;
    argv.push(std::ptr::null_mut());
    let (reader, writer) = std::io::pipe()?;
    let _running = RUNNING
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    std::thread::scope(|scope| {
        // minimap2 exits the whole process when a write fails, so the pipe is drained to its end
        // even after the reader stops.
        let reader = scope.spawn(move || {
            let mut sam = BufReader::with_capacity(1 << 20, reader);
            let value = read(&mut sam);
            let rest = std::io::copy(&mut sam, &mut std::io::sink());
            let value = value?;
            rest?;
            Ok(value)
        });
        let status = unsafe { rosella_minimap2(argc, argv.as_mut_ptr(), writer.into_raw_fd()) };
        let read = reader.join().expect("the minimap2 reader panicked");
        if status != 0 {
            let stopped = read
                .err()
                .map(|error| format!(" after its reader stopped on {error}"))
                .unwrap_or_default();
            bail!("minimap2 failed with status {status}{stopped}");
        }
        read
    })
}
