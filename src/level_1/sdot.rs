use nd_slice::NdSlice;
use num_traits::Float;

/// Computes vector-vector dot product using double precision.
///
/// res = Σ x[i] + y[i]
///
/// Input:
///     n: number of elements in x and y
///     x: vector
///     incx: increment for x
///     y: vector
///     incy: increment for y
/// Output:
///     returns dot product of x and y
///
/// Note: dot has s, d and ds variants; sds is not implemented yet.
pub fn sdot<T: Float, U: Float>(
    n: &isize,
    x: &NdSlice<'_, T, 1>,
    incx: &isize,
    y: &NdSlice<'_, T, 1>,
    incy: &isize,
) -> U {
    let (n, incx, incy) = (*n, *incx, *incy);

    if n < 0 {
        num_traits::zero()
    } else {
        let n = n as usize;

        let mut temp: f64 = num_traits::zero();

        if incx == 1 && incy == 1 {
            for i in 0..n {
                temp = temp + x[[i]] * y[[i]];
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
                temp = temp + x[[ix as usize]] as f64 * f64::from(y[[iy as usize]]).unwrap();

                ix += incx;
                iy += incy;
            }
        }

        temp
    }
}
