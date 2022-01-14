use nd_slice::NdSlice;
use num_traits::Float;

/// Computes vector-vector dot product
///
/// res = Σ x[i] * y[i]
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
pub fn dot<T: Float, U: Float>(
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

        let mut temp = num_traits::zero();

        if incx == 1 && incy == 1 {
            for i in 0..n {
                temp = temp + U::from(x[[i]]).unwrap() * U::from(y[[i]]).unwrap();
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
                temp =
                    temp + U::from(x[[ix as usize]]).unwrap() * U::from(y[[iy as usize]]).unwrap();

                ix += incx;
                iy += incy;
            }
        }

        temp
    }
}

/// Computes vector-vector dot product with double precision.
/// The two vectors are cast to double precision floats and a constant sb is added to the result.
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
pub fn sdsdot(
    n: &isize,
    sb: &f32,
    x: &NdSlice<'_, f32, 1>,
    incx: &isize,
    y: &NdSlice<'_, f32, 1>,
    incy: &isize,
) -> f32 {
    let (n, incx, incy) = (*n, *incx, *incy);

    if n < 0 {
        *sb
    } else {
        let n = n as usize;

        let mut temp = *sb;

        if incx == 1 && incy == 1 {
            for i in 0..n {
                temp = temp + (x[[i]] as f64 * y[[i]] as f64) as f32;
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
                temp = temp + (x[[ix as usize]] as f64 * y[[iy as usize]] as f64) as f32;

                ix += incx;
                iy += incy;
            }
        }

        temp
    }
}
