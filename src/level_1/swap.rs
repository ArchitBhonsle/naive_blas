use nd_slice::NdSliceMut;

/// Given two vectors x and y, swap their elements.
///
/// Input:
///     n: number of elements in x and y
///     x: array
///     incx: increment for x
///     y: array
///     incy: increment for x
/// Output:
///     x: modified x
///     y: modified y
pub fn swap<T>(
    n: isize,
    mut x: NdSliceMut<'_, T, 1>,
    incx: isize,
    mut y: NdSliceMut<'_, T, 1>,
    incy: isize,
) {
    if n < 0 {
        return;
    }
    let n = n as usize;

    if incx == 1 && incy == 1 {
        for i in 0..n {
            std::mem::swap(&mut x[[i]], &mut y[[i]]);
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
            std::mem::swap(&mut x[[(ix as usize)]], &mut y[[(iy as usize)]]);

            ix += incx;
            iy += incy;
        }
    }
}
