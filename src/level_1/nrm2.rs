use nd_slice::NdSlice;
use num_complex::Complex;
use num_traits::Float;

/// Computes the Euclidean norm of a vector.
///
/// res = ‖ x ‖
///
/// Input:
///     n: number of elements in x and y
///     x: vector
///     incx: increment for x
/// Output:
///     euclidean norm of the vector
pub fn nrm2<F: Float>(n: &isize, x: &mut NdSlice<'_, F, 1>, incx: &isize) -> F {
    let (n, incx) = (*n, *incx);

    let (zero, one): (F, F) = (num_traits::zero(), num_traits::one());

    if n < 1 || incx < 1 {
        zero
    } else if n == 1 {
        x[[0]]
    } else {
        let (n, incx) = (n as usize, incx as usize);
        let (mut scale, mut ssq) = (zero, one);

        for ix in (0..(n - 1) * incx).step_by(incx) {
            if x[[ix]] != zero {
                let absxi = x[[ix]].abs();

                if scale < absxi {
                    ssq = one + ssq * (scale / absxi).powi(2);
                    scale = absxi;
                } else {
                    ssq = ssq + (absxi / scale).powi(2);
                }
            }
        }

        scale * ssq.sqrt()
    }
}

pub fn real_complex_nrm2<F: Float>(n: &isize, x: NdSlice<'_, Complex<F>, 1>, incx: &isize) -> F {
    let (n, incx) = (*n, *incx);

    let (zero, one): (F, F) = (num_traits::zero(), num_traits::one());

    if n < 1 || incx < 1 {
        zero
    } else {
        let (n, incx) = (n as usize, incx as usize);
        let (mut scale, mut ssq) = (zero, one);

        for ix in (0..(n - 1) * incx).step_by(incx) {
            if x[[ix]].re != zero {
                let abs_xi_re = x[[ix]].re.abs();

                if scale < abs_xi_re {
                    ssq = one + ssq * (scale / abs_xi_re).powi(2);
                    scale = abs_xi_re;
                } else {
                    ssq = ssq + (abs_xi_re / scale).powi(2);
                }
            }

            if x[[ix]].im != zero {
                let abs_xi_im = x[[ix]].im.abs();

                if scale < abs_xi_im {
                    ssq = one + ssq * (scale / abs_xi_im).powi(2);
                    scale = abs_xi_im;
                } else {
                    ssq = ssq + (abs_xi_im / scale).powi(2);
                }
            }
        }

        scale * ssq.sqrt()
    }
}
