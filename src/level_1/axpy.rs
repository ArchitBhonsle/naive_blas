use nd_slice::NdSliceMut;
use num_complex::Complex;
use num_traits::Float;

/// Computes a vector-scalar product and adds the result to the vector.
///
/// y = ax + y
///
/// Input:
///     n: number of elements in x and y
///     a: scalar to multiply with
///     x: vector
///     incx: increment for x
///     y: vector
///     incy: increment for y
/// Output:
///     y: modified y
pub fn axpy_real<T: Float>(
    n: &isize,
    a: &T,
    x: &mut NdSliceMut<'_, T, 1>,
    incx: &isize,
    y: &mut NdSliceMut<'_, T, 1>,
    incy: &isize,
) {
    let (n, a, incx, incy) = (*n, *a, *incx, *incy);

    if n < 0 {
        return;
    }
    let n = n as usize;

    if a == num_traits::zero() {
        return;
    }

    if incx == 1 && incy == 1 {
        for i in 0..n {
            y[[i]] = y[[i]] + a * x[[i]];
        }
    } else {
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
            y[[iy as usize]] = y[[iy as usize]] + a * x[[ix as usize]];

            ix += incx;
            iy += incy;
        }
    }
}

/// Computes a vector-scalar product and adds the result to the vector.
///
/// y = ax + y
///
/// Input:
///     n: number of elements in x and y
///     a: scalar to multiply with
///     x: vector
///     incx: increment for x
///     y: vector
///     incy: increment for y
/// Output:
///     y: modified y
pub fn axpy_complex<T: Float>(
    n: &isize,
    a: &Complex<T>,
    x: &mut NdSliceMut<'_, Complex<T>, 1>,
    incx: &isize,
    y: &mut NdSliceMut<'_, Complex<T>, 1>,
    incy: &isize,
) {
    let (n, incx, incy) = (*n, *incx, *incy);
    if n < 0 {
        return;
    }
    let n = n as usize;

    if a.re.abs() + a.im.abs() == num_traits::zero() {
        return;
    }

    if incx == 1 && incy == 1 {
        for i in 0..n {
            y[[i]] = y[[i]] + a * x[[i]];
        }
    } else {
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
            y[[iy as usize]] = y[[iy as usize]] + a * x[[ix as usize]];

            ix += incx;
            iy += incy;
        }
    }
}
