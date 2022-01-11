use nd_slice::NdSliceMut;
use num_traits::{Float, NumCast};

/// Performs rotation of points in a plane.
///
/// Given two vectors x and y, each element is replaced as:
///     x[i] = c*x[i] + s*y[i]
///     y[i] = c*y[i] + s*x[i]
///
/// Input:
///     n: number of elements in each of the vectors
///     x: array
///     incx: increment for x
///     y: array
///     incy: increment for x
///     c: scalar
///     s: scalar
/// Output:
///     x: modified x
///     y: modified y
///
/// Note: Intel MKL also defines rotg for c and z too.
pub fn rot<T: Float>(
    n: &isize,
    x: &mut NdSliceMut<'_, T, 1>,
    incx: &isize,
    y: &mut NdSliceMut<'_, T, 1>,
    incy: &isize,
    c: &isize,
    s: &isize,
) {
    let (n, incx, incy) = (*n, *incx, *incy);
    let (c, s): (T, T) = (NumCast::from(*c).unwrap(), NumCast::from(*s).unwrap());

    if n < 0 {
        return;
    }
    let n = n as usize;

    if incx == 1 && incy == 1 {
        for i in 0..n {
            let temp = c * x[[i]] + s * y[[i]];
            y[[i]] = c * y[[i]] - s * x[[i]];
            x[[i]] = temp;
        }

        return;
    }

    let mut ix = if incx < 0 {
        (1 - (n as isize)) * incx
    } else {
        0
    };
    let mut iy = if incy < 0 {
        (1 - (n as isize)) * incy
    } else {
        0
    };

    for _ in 0..n {
        let (ix_inner, iy_inner) = (ix as usize, iy as usize);

        let temp = c * x[[ix_inner]] + s * y[[iy_inner]];
        y[[iy_inner]] = c * y[[iy_inner]] - s * x[[ix_inner]];
        x[[ix_inner]] = temp;

        ix += incx;
        iy += incy;
    }
}
