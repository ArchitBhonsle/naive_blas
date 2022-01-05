use nd_slice::NdSliceMut;
use num_traits::{Float, NumCast};

/// Performs rotation of points in a plane
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
///     x: rotated x
///     y: rotated y
pub fn rot<T: Float>(
    n: isize,
    mut x: NdSliceMut<'_, T, 1>,
    incx: isize,
    mut y: NdSliceMut<'_, T, 1>,
    incy: isize,
    c: isize,
    s: isize,
) {
    if n < 0 {
        return;
    }
    let n = n as usize;
    let c: T = NumCast::from(c).unwrap();
    let s: T = NumCast::from(s).unwrap();

    if incx == 1 && incy == 1 {
        for i in 0..n {
            let temp = c * x[[i]] + s * y[[i]];
            y[[i]] = c * y[[i]] - s * x[[i]];
            x[[i]] = temp;
        }

        return;
    }

    let mut ix = if incx < 0 {
        (-(n as isize) + 1) * incx + 1
    } else {
        1
    };
    let mut iy = if incy < 0 {
        (-(n as isize) + 1) * incy + 1
    } else {
        1
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
