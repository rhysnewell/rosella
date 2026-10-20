use log::{info, warn};
use ndarray::{Array2, ArrayView2, ArrayViewMut2};

use crate::learned::device::{self, Found};

pub struct Gpu;

// Called on the thread that trains, which holds the device until the guard drops.
pub fn gpu(wanted: bool) -> Option<Gpu> {
    if !wanted {
        return None;
    }
    let opened = match device::open() {
        Ok(Found::Device(name)) => Ok(name),
        Ok(Found::Missing(why)) => Err(why),
        Err(why) => Err(format!("{why:#}")),
    };
    match opened {
        Ok(name) => {
            info!("GPU: {name}");
            Some(Gpu)
        }
        Err(why) => {
            warn!("GPU: {why}, training on the CPU");
            None
        }
    }
}

impl Drop for Gpu {
    fn drop(&mut self) {
        device::close();
    }
}

struct Out {
    at: *mut f32,
    dim: (usize, usize),
    strides: [isize; 2],
}

pub fn product(
    alpha: f32,
    a: ArrayView2<f32>,
    b: ArrayView2<f32>,
    beta: f32,
    mut c: ArrayViewMut2<f32>,
) {
    let out = Out {
        at: c.as_mut_ptr(),
        dim: c.dim(),
        strides: [c.strides()[0], c.strides()[1]],
    };
    unsafe { start(alpha, &a, &b, beta, out) };
}

// With beta zero the output is never read, so a fresh table need not be zeroed first.
pub fn fresh(alpha: f32, a: ArrayView2<f32>, b: ArrayView2<f32>) -> Array2<f32> {
    let mut table = Array2::<f32>::uninit((a.nrows(), b.ncols()));
    let dim = table.dim();
    let out = Out {
        at: table.as_mut_ptr().cast(),
        dim,
        strides: [dim.1 as isize, 1],
    };
    unsafe { start(alpha, &a, &b, 0.0, out) };
    unsafe { table.assume_init() }
}

#[cfg(target_os = "macos")]
#[derive(Clone, Copy)]
struct Layout {
    stride: usize,
    turned: bool,
}

#[cfg(target_os = "macos")]
impl Layout {
    fn of((rows, columns): (usize, usize), strides: &[isize]) -> Option<Self> {
        let fits = |step: isize, along: usize, span: usize| {
            along <= 1 || usize::try_from(step).is_ok_and(|step| step >= span)
        };
        let plain = (columns <= 1 || strides[1] == 1) && fits(strides[0], rows, columns);
        let turned = (rows <= 1 || strides[0] == 1) && fits(strides[1], columns, rows);
        match (plain, turned) {
            (true, _) => Some(Self {
                stride: if rows > 1 {
                    strides[0] as usize
                } else {
                    columns
                },
                turned: false,
            }),
            (false, true) => Some(Self {
                stride: if columns > 1 {
                    strides[1] as usize
                } else {
                    rows
                },
                turned: true,
            }),
            _ => None,
        }
    }
}

#[cfg(target_os = "macos")]
unsafe fn start(alpha: f32, a: &ArrayView2<f32>, b: &ArrayView2<f32>, beta: f32, out: Out) {
    use cblas_sys::{CblasNoTrans, CblasRowMajor, CblasTrans, cblas_sgemm};
    let ((m, n), k) = (out.dim, a.ncols());
    if m == 0 {
        return;
    }
    let layout = |dim, strides: &[isize]| {
        Layout::of(dim, strides).expect("products take plain or turned tables")
    };
    let (left, right) = (layout(a.dim(), a.strides()), layout(b.dim(), b.strides()));
    let target = layout(out.dim, &out.strides);
    assert!(!target.turned, "products write plain tables");
    let turn = |layout: Layout| match layout.turned {
        true => CblasTrans,
        false => CblasNoTrans,
    };
    unsafe {
        cblas_sgemm(
            CblasRowMajor,
            turn(left),
            turn(right),
            m as i32,
            n as i32,
            k as i32,
            alpha,
            a.as_ptr(),
            left.stride as i32,
            b.as_ptr(),
            right.stride as i32,
            beta,
            out.at,
            target.stride as i32,
        );
    }
}

// The gemm crate sums each output cell in one order, so the bits do not depend on the thread
// count. It ran the training shapes 1.6 times faster than matrixmultiply on 4 AVX2 cores.
#[cfg(not(target_os = "macos"))]
unsafe fn start(alpha: f32, a: &ArrayView2<f32>, b: &ArrayView2<f32>, beta: f32, out: Out) {
    let ((m, n), k) = (out.dim, a.ncols());
    let ([ra, ca], [rb, cb], [rc, cc]) = (
        [a.strides()[0], a.strides()[1]],
        [b.strides()[0], b.strides()[1]],
        out.strides,
    );
    unsafe {
        gemm::gemm(
            m,
            n,
            k,
            out.at,
            cc,
            rc,
            beta != 0.0,
            a.as_ptr(),
            ca,
            ra,
            b.as_ptr(),
            cb,
            rb,
            beta,
            alpha,
            false,
            false,
            false,
            gemm::Parallelism::Rayon(rayon::current_num_threads()),
        );
    }
}
