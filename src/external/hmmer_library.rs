use std::ffi::{CStr, CString, c_char, c_int};
use std::os::unix::ffi::OsStrExt;
use std::path::Path;
use std::ptr::NonNull;

use anyhow::{Result, anyhow, bail};

#[repr(C)]
struct RawModels {
    _held: [u8; 0],
}

#[repr(C)]
struct RawTargets {
    _held: [u8; 0],
}

unsafe extern "C" {
    fn rosella_models_read(path: *const c_char, error: *mut c_char, room: usize) -> *mut RawModels;
    fn rosella_models_count(models: *const RawModels) -> c_int;
    fn rosella_models_free(models: *mut RawModels);
    fn rosella_targets_new() -> *mut RawTargets;
    fn rosella_targets_add(
        targets: *mut RawTargets,
        models: *const RawModels,
        id: u64,
        protein: *const c_char,
        length: i64,
    ) -> c_int;
    fn rosella_targets_residues(targets: *const RawTargets) -> i64;
    fn rosella_targets_clear(targets: *mut RawTargets);
    fn rosella_targets_free(targets: *mut RawTargets);
    fn rosella_search(
        models: *const RawModels,
        model: c_int,
        targets: *const RawTargets,
        floor: f64,
        size: *mut usize,
    ) -> *mut c_char;
    fn rosella_text_free(text: *mut c_char);
}

pub struct Models(NonNull<RawModels>);

pub struct Targets(NonNull<RawTargets>);

// A search only reads the models and the held proteins, so any number of threads may share them.
unsafe impl Send for Models {}
unsafe impl Sync for Models {}
unsafe impl Send for Targets {}
unsafe impl Sync for Targets {}

impl Models {
    pub fn read(path: &Path) -> Result<Self> {
        let name = CString::new(path.as_os_str().as_bytes())?;
        let mut error = [0 as c_char; 256];
        let raw = unsafe { rosella_models_read(name.as_ptr(), error.as_mut_ptr(), error.len()) };
        NonNull::new(raw).map(Self).ok_or_else(|| {
            let reason = unsafe { CStr::from_ptr(error.as_ptr()) }.to_string_lossy();
            anyhow!("HMMER could not read {}: {reason}", path.display())
        })
    }

    pub fn len(&self) -> usize {
        unsafe { rosella_models_count(self.0.as_ptr()) as usize }
    }

    pub fn search(&self, model: usize, targets: &Targets, floor: f64) -> Result<String> {
        let mut size = 0;
        let text = unsafe {
            rosella_search(
                self.0.as_ptr(),
                model as c_int,
                targets.0.as_ptr(),
                floor,
                &mut size,
            )
        };
        let text = NonNull::new(text).ok_or_else(|| anyhow!("HMMER failed on model {model}"))?;
        let bytes =
            unsafe { std::slice::from_raw_parts(text.as_ptr().cast::<u8>(), size) }.to_vec();
        unsafe { rosella_text_free(text.as_ptr()) };
        Ok(String::from_utf8(bytes)?)
    }
}

impl Drop for Models {
    fn drop(&mut self) {
        unsafe { rosella_models_free(self.0.as_ptr()) }
    }
}

impl Targets {
    pub fn new() -> Result<Self> {
        let raw = unsafe { rosella_targets_new() };
        NonNull::new(raw)
            .map(Self)
            .ok_or_else(|| anyhow!("HMMER could not hold proteins"))
    }

    pub fn add(&mut self, models: &Models, id: usize, protein: &str) -> Result<()> {
        let status = unsafe {
            rosella_targets_add(
                self.0.as_ptr(),
                models.0.as_ptr(),
                id as u64,
                protein.as_ptr().cast::<c_char>(),
                protein.len() as i64,
            )
        };
        if status != 0 {
            bail!("HMMER refused protein {id}, status {status}");
        }
        Ok(())
    }

    pub fn residues(&self) -> usize {
        unsafe { rosella_targets_residues(self.0.as_ptr()) as usize }
    }

    pub fn clear(&mut self) {
        unsafe { rosella_targets_clear(self.0.as_ptr()) }
    }
}

impl Drop for Targets {
    fn drop(&mut self) {
        unsafe { rosella_targets_free(self.0.as_ptr()) }
    }
}
