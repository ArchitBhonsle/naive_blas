use nd_slice::NdSliceMut;
use num_traits::Float;

/// Performs modified Givens rotation of points in the plane.
///
/// Given two vectors x and y, replace their elements as:
///     | x[i] | = H | x[i] |
///     | y[i] |     | y[i] |
///
/// Input:
///     n: number of elements in each of the vectors
///     x: array
///     incx: increment for x
///     y: array
///     incy: increment for x
///     param: array of size 5
///         param[0] -> flag, param[1-4] -> h11, h21, h12, h22
///         if flag = -1 -> H = |  h11  h12 |
///                             |  h21  h22 |
///         if flag =  0 -> H = |  1.0  h12 |
///                             |  h21  1.0 |
///         if flag = -1 -> H = |  h11  1.0 |
///                             | -1.0  h22 |
///         if flag = -1 -> H = |  1.0  0.0 |
///                             |  0.0  1.0 |
///         Note: -1.0, 0.0, 1.0 are assumed based on the flag and may not be set in param
/// Output:
///     x: modified x
///     y: modified y
pub fn rotm<T: Float>(
    n: isize,
    mut x: NdSliceMut<'_, T, 1>,
    incx: isize,
    mut y: NdSliceMut<'_, T, 1>,
    incy: isize,
    param: NdSliceMut<'_, T, 1>,
) {
    let zero: T = num_traits::zero();
    let two: T = num_traits::one::<T>() + num_traits::one();

    let flag = param[[0]];

    if n < 0 || flag + two == zero {
        return;
    }
    let n = n as usize;

    if incx == incy && incx > 0 {
        let incx = incx as usize;
        let nsteps = n * incx;

        if flag < zero {
            let (h11, h12, h21, h22) = (param[[1]], param[[3]], param[[2]], param[[4]]);

            for i in (0..nsteps).step_by(incx) {
                let w = x[[i]];
                let z = y[[i]];

                x[[i]] = w * h11 + z * h12;
                y[[i]] = w * h21 + z * h22;
            }
        } else if flag == zero {
            let (h12, h21) = (param[[3]], param[[2]]);

            for i in (0..nsteps).step_by(incx) {
                let w = x[[i]];
                let z = y[[i]];

                x[[i]] = w + z * h12;
                y[[i]] = w * h21 + z;
            }
        } else {
            let (h11, h22) = (param[[1]], param[[4]]);

            for i in (0..nsteps).step_by(incx) {
                let w = x[[i]];
                let z = y[[i]];

                x[[i]] = w * h11 + z;
                y[[i]] = -w + z * h22;
            }
        }
    } else {
        let mut kx = if incx < 0 {
            (1 - (n as isize)) * incx
        } else {
            0
        };
        let mut ky = if incy < 0 {
            (1 - (n as isize)) * incy
        } else {
            0
        };

        if flag < zero {
            let (h11, h12, h21, h22) = (param[[1]], param[[3]], param[[2]], param[[4]]);

            for _ in 0..n {
                let (kx_inner, ky_inner) = (kx as usize, ky as usize);

                let w = x[[kx_inner]];
                let z = y[[ky_inner]];

                x[[kx_inner]] = w * h11 + z * h12;
                y[[ky_inner]] = w * h21 + z * h22;

                kx += incx;
                ky += incy;
            }
        } else if flag == zero {
            let (h12, h21) = (param[[3]], param[[2]]);

            for _ in 0..n {
                let (kx_inner, ky_inner) = (kx as usize, ky as usize);

                let w = x[[kx_inner]];
                let z = y[[ky_inner]];

                x[[kx_inner]] = w + z * h12;
                y[[ky_inner]] = w * h21 + z;

                kx += incx;
                ky += incy;
            }
        } else {
            let (h11, h22) = (param[[1]], param[[4]]);

            for _ in 0..n {
                let (kx_inner, ky_inner) = (kx as usize, ky as usize);

                let w = x[[kx_inner]];
                let z = y[[ky_inner]];

                x[[kx_inner]] = w * h11 + z;
                y[[ky_inner]] = -w + z * h22;

                kx += incx;
                ky += incy;
            }
        }
    }
}
