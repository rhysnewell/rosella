use std::marker::PhantomData;
use std::mem::MaybeUninit;

#[cfg(target_os = "macos")]
use log::info;
use log::warn;
use ndarray::{Array2, ArrayView2, ArrayViewMut2};

#[cfg(target_os = "macos")]
use crate::learned::metal;

pub struct Gpu;

// Called on the thread that trains, which holds the device until the guard drops.
pub fn gpu(wanted: bool) -> Option<Gpu> {
    if !wanted {
        return None;
    }
    #[cfg(target_os = "macos")]
    if let Some(name) = metal::open() {
        info!("GPU: {name}");
        return Some(Gpu);
    }
    warn!("GPU: none found");
    None
}

impl Drop for Gpu {
    fn drop(&mut self) {
        #[cfg(target_os = "macos")]
        metal::close();
    }
}

struct Out {
    at: *mut f32,
    dim: (usize, usize),
    strides: [isize; 2],
}

impl Out {
    fn fresh(table: &mut Array2<MaybeUninit<f32>>) -> Self {
        let dim = table.dim();
        Self {
            at: table.as_mut_ptr().cast(),
            dim,
            strides: [dim.1 as isize, 1],
        }
    }
}

#[cfg(target_os = "macos")]
type Running = metal::Running;
#[cfg(not(target_os = "macos"))]
type Running = std::convert::Infallible;

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
    drop(unsafe { start(alpha, &a, &b, beta, out, false) });
}

// With beta zero the output is never read, so a fresh table need not be zeroed first.
pub fn fresh(alpha: f32, a: ArrayView2<f32>, b: ArrayView2<f32>) -> Array2<f32> {
    let mut table = Array2::uninit((a.nrows(), b.ncols()));
    drop(unsafe { start(alpha, &a, &b, 0.0, Out::fresh(&mut table), false) });
    unsafe { table.assume_init() }
}

// The work is dropped before the table, so the table outlives any GPU still writing it.
pub struct Later<'a> {
    running: Option<Running>,
    table: Array2<MaybeUninit<f32>>,
    held: Option<Array2<f32>>,
    inputs: PhantomData<ArrayView2<'a, f32>>,
}

impl Later<'_> {
    pub fn wait(self) -> Array2<f32> {
        let Later {
            running,
            table,
            held,
            ..
        } = self;
        drop(running);
        drop(held);
        unsafe { table.assume_init() }
    }
}

pub fn later<'a>(alpha: f32, a: ArrayView2<'a, f32>, b: ArrayView2<'a, f32>) -> Later<'a> {
    let mut table = Array2::uninit((a.nrows(), b.ncols()));
    let running = unsafe { start(alpha, &a, &b, 0.0, Out::fresh(&mut table), true) };
    Later {
        running,
        table,
        held: None,
        inputs: PhantomData,
    }
}

// Takes the left side, turned, so a gradient can be handed off while the GPU still reads it.
pub fn later_turned<'a>(alpha: f32, a: Array2<f32>, b: ArrayView2<'a, f32>) -> Later<'a> {
    let mut table = Array2::uninit((a.ncols(), b.ncols()));
    let running = unsafe { start(alpha, &a.t(), &b, 0.0, Out::fresh(&mut table), true) };
    Later {
        running,
        table,
        held: Some(a),
        inputs: PhantomData,
    }
}

#[cfg(target_os = "macos")]
#[derive(Clone, Copy)]
pub(crate) struct Layout {
    pub(crate) rows: usize,
    pub(crate) columns: usize,
    pub(crate) stride: usize,
    pub(crate) turned: bool,
}

#[cfg(target_os = "macos")]
impl Layout {
    fn of(dim: (usize, usize), strides: [isize; 2]) -> Option<Self> {
        let (rows, columns) = dim;
        let fits = |step: isize, along: usize, span: usize| {
            along <= 1 || usize::try_from(step).is_ok_and(|step| step >= span)
        };
        let plain = (columns <= 1 || strides[1] == 1) && fits(strides[0], rows, columns);
        let turned = (rows <= 1 || strides[0] == 1) && fits(strides[1], columns, rows);
        match (plain, turned) {
            (true, _) => Some(Self {
                rows,
                columns,
                stride: if rows > 1 {
                    strides[0] as usize
                } else {
                    columns
                },
                turned: false,
            }),
            (false, true) => Some(Self {
                rows: columns,
                columns: rows,
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

    pub(crate) fn elements(&self) -> usize {
        (self.rows - 1) * self.stride + self.columns
    }

    pub(crate) fn first(self, rows: usize) -> Self {
        match self.turned {
            false => Self { rows, ..self },
            true => Self {
                columns: rows,
                ..self
            },
        }
    }

    fn past(self, rows: usize) -> (usize, Self) {
        let left = match self.turned {
            false => self.rows - rows,
            true => self.columns - rows,
        };
        let offset = if self.turned {
            rows
        } else {
            rows * self.stride
        };
        (offset, self.first(left))
    }
}

#[cfg(target_os = "macos")]
type Ends = (
    (*const f32, Layout),
    (*const f32, Layout),
    (*mut f32, Layout),
);

#[cfg(target_os = "macos")]
unsafe fn start(
    alpha: f32,
    a: &ArrayView2<f32>,
    b: &ArrayView2<f32>,
    beta: f32,
    out: Out,
    any_size: bool,
) -> Option<Running> {
    let layout = |dim, strides: &[isize]| {
        Layout::of(dim, [strides[0], strides[1]]).expect("products take plain or turned tables")
    };
    let (left, right) = (layout(a.dim(), a.strides()), layout(b.dim(), b.strides()));
    let target = layout(out.dim, &out.strides);
    assert!(!target.turned, "products write plain tables");
    let shape = (out.dim.0, out.dim.1, a.ncols());
    let ends = ((a.as_ptr(), left), (b.as_ptr(), right), (out.at, target));
    if any_size {
        if let Some(running) = metal::start(alpha, ends, beta, shape) {
            return Some(running);
        }
    } else if let Some(split) = metal::split(alpha, ends, beta, shape) {
        let rows = split.rows();
        let ((a_offset, left), (c_offset, target)) = (left.past(rows), target.past(rows));
        let rest = (
            (unsafe { a.as_ptr().add(a_offset) }, left),
            (b.as_ptr(), right),
            (unsafe { out.at.add(c_offset) }, target),
        );
        unsafe { matrix_unit(alpha, rest, beta, (shape.0 - rows, shape.1, shape.2)) };
        split.finish();
        return None;
    }
    unsafe { matrix_unit(alpha, ends, beta, shape) };
    None
}

// The same call ndarray makes into Accelerate, so the sums come out the same.
#[cfg(target_os = "macos")]
unsafe fn matrix_unit(alpha: f32, ends: Ends, beta: f32, (m, n, k): (usize, usize, usize)) {
    use cblas_sys::{CBLAS_TRANSPOSE, CblasNoTrans, CblasRowMajor, CblasTrans, cblas_sgemm};
    if m == 0 {
        return;
    }
    let turn = |layout: Layout| -> CBLAS_TRANSPOSE {
        if layout.turned {
            CblasTrans
        } else {
            CblasNoTrans
        }
    };
    let ((a, left), (b, right), (c, target)) = ends;
    unsafe {
        cblas_sgemm(
            CblasRowMajor,
            turn(left),
            turn(right),
            m as i32,
            n as i32,
            k as i32,
            alpha,
            a,
            left.stride as i32,
            b,
            right.stride as i32,
            beta,
            c,
            target.stride as i32,
        );
    }
}

// One product runs on one core, so the longer side of the output is shared out. Each cell is
// still summed by one call in one order, so the result does not depend on the thread count.
#[cfg(not(target_os = "macos"))]
unsafe fn start(
    alpha: f32,
    a: &ArrayView2<f32>,
    b: &ArrayView2<f32>,
    beta: f32,
    out: Out,
    _: bool,
) -> Option<Running> {
    use rayon::prelude::*;
    let ((m, k), n) = (a.dim(), b.ncols());
    let ([ra, ca], [rb, cb], [rc, cc]) = (
        [a.strides()[0], a.strides()[1]],
        [b.strides()[0], b.strides()[1]],
        out.strides,
    );
    let (a_at, b_at, c_at) = (a.as_ptr() as usize, b.as_ptr() as usize, out.at as usize);
    let pieces = 2 * rayon::current_num_threads();
    let by_rows = m >= n;
    let length = if by_rows { m } else { n };
    let step = length.div_ceil(pieces).max(1);
    (0..length.div_ceil(step))
        .into_par_iter()
        .for_each(|piece| {
            let start = piece * step;
            let span = step.min(length - start);
            let (rows, columns) = if by_rows { (span, n) } else { (m, span) };
            let (a_from, b_from, c_from) = if by_rows {
                (start as isize * ra, 0, start as isize * rc)
            } else {
                (0, start as isize * cb, start as isize * cc)
            };
            unsafe {
                matrixmultiply::sgemm(
                    rows,
                    k,
                    columns,
                    alpha,
                    (a_at as *const f32).offset(a_from),
                    ra,
                    ca,
                    (b_at as *const f32).offset(b_from),
                    rb,
                    cb,
                    beta,
                    (c_at as *mut f32).offset(c_from),
                    rc,
                    cc,
                );
            }
        });
    None
}
