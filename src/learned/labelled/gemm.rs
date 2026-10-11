use ndarray::{ArrayView2, ArrayViewMut2, Axis};
use rayon::prelude::*;

// One product runs on one core, so the longer side of the output is shared out. Each cell is
// still summed by one call in one order, so the result does not depend on the thread count.
pub fn product(
    alpha: f32,
    a: ArrayView2<f32>,
    b: ArrayView2<f32>,
    beta: f32,
    mut c: ArrayViewMut2<f32>,
) {
    let (rows, columns) = c.dim();
    let pieces = 2 * rayon::current_num_threads();
    if rows >= columns {
        let step = rows.div_ceil(pieces).max(1);
        c.axis_chunks_iter_mut(Axis(0), step)
            .zip(a.axis_chunks_iter(Axis(0), step))
            .collect::<Vec<_>>()
            .into_par_iter()
            .for_each(|(mut c, a)| cores(alpha, &a, &b, beta, &mut c));
    } else {
        let step = columns.div_ceil(pieces).max(1);
        c.axis_chunks_iter_mut(Axis(1), step)
            .zip(b.axis_chunks_iter(Axis(1), step))
            .collect::<Vec<_>>()
            .into_par_iter()
            .for_each(|(mut c, b)| cores(alpha, &a, &b, beta, &mut c));
    }
}

// What ndarray ran before it was built with BLAS. On macOS it now hands products to Accelerate.
fn cores(
    alpha: f32,
    a: &ArrayView2<f32>,
    b: &ArrayView2<f32>,
    beta: f32,
    c: &mut ArrayViewMut2<f32>,
) {
    let ((m, k), n) = (a.dim(), b.ncols());
    assert_eq!((b.nrows(), c.dim()), (k, (m, n)), "product shapes agree");
    unsafe {
        matrixmultiply::sgemm(
            m,
            k,
            n,
            alpha,
            a.as_ptr(),
            a.strides()[0],
            a.strides()[1],
            b.as_ptr(),
            b.strides()[0],
            b.strides()[1],
            beta,
            c.as_mut_ptr(),
            c.strides()[0],
            c.strides()[1],
        );
    }
}
