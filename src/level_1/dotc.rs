use nd_slice::NdSlice;
use num_complex::Complex;
use num_traits::Float;

/// Computes dot product of conjugated vector with another vector
///
/// res = Σ conjugate(x[i]) * y[i]
///
/// Input:
///     n: number of elements in x and y
///     x: vector
///     incx: increment for x
///     y: vector
///     incy: increment for y
/// Output:
///     returns dot product of complex vectors x and y
pub fn dotc<T: Float>(
    n: &isize,
    x: &mut NdSlice<'_, Complex<T>, 1>,
    incx: &isize,
    y: &mut NdSlice<'_, Complex<T>, 1>,
    incy: &isize,
) -> Complex<T> {
    let (n, incx, incy) = (*n, *incx, *incy);

    if n < 0 {
        Complex::new(num_traits::zero(), num_traits::zero())
    } else {
        let n = n as usize;

        let mut temp = Complex::new(num_traits::zero(), num_traits::zero());

        if incx == 1 && incy == 1 {
            for i in 0..n {
                temp = temp + x[[i]].conj() * y[[i]];
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
                temp = temp + x[[ix as usize]].conj() * y[[iy as usize]];

                ix += incx;
                iy += incy;
            }
        }

        temp
    }
}
