use nd_slice::NdSliceMut;

/// Copies one vector to another.
///
/// y = x
///
/// Input:
///     n: number of elements in x and y
///     x: array
///     incx: increment for x
///     y: array
///     incy: increment for x
/// Output:
///     y: copies x
pub fn copy<T: Copy /* both Float and Complex<T> implement Copy */>(
    n: &isize,
    x: &mut NdSliceMut<'_, T, 1>,
    incx: &isize,
    y: &mut NdSliceMut<'_, T, 1>,
    incy: &isize,
) {
    let (n, incx, incy) = (*n, *incx, *incy);
    if n < 0 {
        return;
    }
    let n = n as usize;

    if incx == 1 && incy == 1 {
        for i in 0..n {
            y[[i]] = x[[i]];
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
            y[[ix as usize]] = x[[iy as usize]];

            ix += incx;
            iy += incy;
        }
    }
}
