use nd_slice::NdSliceMut;
use num_complex::Complex;
use num_traits::Float;

/// Computes the product of a vector by a scalar.
///
/// x = a * x
///
/// Input:
///     n: number of elements in the vector
///     a: scalar
///     x: vector
///     incx: increment for x
/// Output:
///     x: modified x
pub fn scal_real<T: Float>(n: &isize, a: &T, x: &mut NdSliceMut<'_, T, 1>, incx: &isize) {
    let (n, a, incx) = (*n, *a, *incx);

    if n < 0 || incx < 0 {
        return;
    }
    let (n, incx) = (n as usize, incx as usize);

    // both branches are actually the same. Netlib does some weird optimization in for incx == 1
    if incx == 1 {
        for i in 0..n {
            x[[i]] = x[[i]] * a;
        }
    } else {
        let nincx = n * incx;
        for i in (0..nincx).step_by(incx) {
            x[[i]] = x[[i]] * a;
        }
    }
}

/// Computes the product of a vector by a scalar
///
/// Input:
///     n: number of elements in the vector
///     a: scalar
///     x: vector
///     incx: increment for x
/// Output:
///     x: modified x
pub fn scal_complex<T: Float>(
    n: isize,
    a: Complex<T>,
    mut x: NdSliceMut<'_, Complex<T>, 1>,
    incx: isize,
) {
    if n < 0 || incx < 0 {
        return;
    }
    let (n, incx) = (n as usize, incx as usize);

    // both branches are actually the same. Netlib does some weird optimization in for incx == 1
    if incx == 1 {
        for i in 0..n {
            x[[i]] = x[[i]] * a;
        }
    } else {
        let nincx = n * incx;
        for i in (0..nincx).step_by(incx) {
            x[[i]] = x[[i]] * a;
        }
    }
}

/// Computes the product of a vector by a scalar
///
/// Input:
///     n: number of elements in the vector
///     a: scalar
///     x: vector
///     incx: increment for x
/// Output:
///     x: modified x
pub fn scal_complex_real<T: Float>(
    n: isize,
    a: T,
    mut x: NdSliceMut<'_, Complex<T>, 1>,
    incx: isize,
) {
    if n < 0 || incx < 0 {
        return;
    }
    let (n, incx) = (n as usize, incx as usize);

    // both branches are actually the same. Netlib does some weird optimization in for incx == 1
    if incx == 1 {
        for i in 0..n {
            x[[i]] = x[[i]].scale(a);
        }
    } else {
        let nincx = n * incx;
        for i in (0..nincx).step_by(incx) {
            x[[i]] = x[[i]].scale(a);
        }
    }
}
