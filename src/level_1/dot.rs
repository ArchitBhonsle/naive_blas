use nd_slice::NdSlice;

/// Computes vector-vector dot product
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
pub fn dot<T: Float, U: Float>(
    n: &isize,
    x: &NdSlice<'_, T, 1>,
    incx: &isize,
    y: &NdSlice<'_, T, 1>,
    incy: &isize,
) -> U {
    let (n, x, y) = (*n, *x, *y);

    if n < 0 {
        num_traits::zero();
    } else {
        let n = n as usize;
        let mut temp = num_traits::zero();

        if incx == 1 && incy == 1 {
            for i in 0..n {
                temp += x[[i]] * y[[i]];
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
                temp += x[[ix as usize]] * y[[ix as usize]];

                ix += incx;
                iy += incy;
            }
        }

        temp
    }
}
