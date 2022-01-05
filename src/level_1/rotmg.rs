use nd_slice::NdSliceMut;
use num_traits::{Float, NumCast};

/// Computes the parameters for a modified Givens rotation.
///
/// Given the co-ordinates (x1, y1), this function calculates the components of a modified Givens
/// transformation matrix H that zeros the y-component of the resulting vector:
///
/// | x1 | = H | x1 √d1 |
/// |  0 | =   | y1 √d2 |
///
/// Input:
///     d1: scaling factor for the x-coordinate of the input vector
///     d2: scaling factor for the y-coordinate of the input vector
///     x1: x-coordinate of the input vector
///     y1: y-coordinate of the input vector
/// Output:
///     d1: first diagonal element of the updated matrix
///     d2: second diagonal element of the updated matrix
///     x1: x-coordinate of the rotated vector before scaling
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
pub fn rotmg<T: Float>(
    d1: &mut T,
    d2: &mut T,
    x1: &mut T,
    y1: &mut T,
    mut param: NdSliceMut<'_, T, 1>,
) {
    let zero: T = num_traits::zero();
    let one: T = num_traits::one();
    let two: T = one + one;

    let gam: T = NumCast::from(4096.).unwrap();
    let gamsq: T = NumCast::from(1.67772e7).unwrap();
    let rgamsq: T = NumCast::from(5.96046e-8).unwrap();

    let mut flag;
    let (mut h11, mut h12, mut h21, mut h22) = (zero, zero, zero, zero);
    if *d1 < zero {
        flag = -one;
        h11 = zero;
        h12 = zero;
        h21 = zero;
        h22 = zero;

        *d1 = zero;
        *d2 = zero;
        *x1 = zero;
    } else {
        let p2 = *d2 * *y1;
        if p2 == zero {
            flag = -two;
            param[[0]] = flag;
            return;
        }

        let p1 = *d1 * *x1;
        let q2 = p2 * *y1;
        let q1 = p1 * *x1;

        if q1.abs() > q2.abs() {
            h21 = -*y1 / *x1;
            h12 = p2 / p1;

            let u = one - h12 * h21;

            if u > zero {
                flag = zero;
                *d1 = *d1 / u;
                *d2 = *d2 / u;
                *x1 = *x1 * u;
            } else {
                // This should only trigger due to rounding errors
                // DOI: 10.1145/355841.355847
                flag = -one;
                h11 = zero;
                h12 = zero;
                h21 = zero;
                h22 = zero;

                *d1 = zero;
                *d2 = zero;
                *x1 = zero;
            }
        } else {
            if q2 < zero {
                flag = -one;
                h11 = zero;
                h12 = zero;
                h21 = zero;
                h22 = zero;

                *d1 = zero;
                *d2 = zero;
                *x1 = zero;
            } else {
                flag = one;
                h11 = p1 / p2;
                h22 = *x1 / *y1;
                let u = one + h11 * h22;
                let temp = *d2 / u;

                *d2 = *d1 / u;
                *d1 = temp;
                *x1 = *y1 * u;
            }
        }

        // scale check
        if *d1 != zero {
            while (*d1 < rgamsq) || (*d1 > gamsq) {
                if flag == zero {
                    h11 = one;
                    h22 = one;
                    flag = -one;
                } else {
                    h21 = -one;
                    h12 = one;
                    flag = -one;
                }

                if *d1 < rgamsq {
                    *d1 = *d1 * gam.powi(2);
                    *x1 = *x1 / gam;
                    h11 = h11 / gam;
                    h12 = h12 / gam;
                } else {
                    *d1 = *d1 / gam.powi(2);
                    *x1 = *x1 * gam;
                    h11 = h11 * gam;
                    h12 = h12 * gam;
                }
            }
        }

        if *d2 != zero {
            while (*d2 < rgamsq) || (*d2 > gamsq) {
                if flag == zero {
                    h11 = one;
                    h22 = one;
                    flag = -one;
                } else {
                    h21 = -one;
                    h12 = one;
                    flag = -one;
                }

                if d2.abs() < rgamsq {
                    *d2 = *d2 * gam.powi(2);
                    h21 = h21 / gam;
                    h22 = h22 / gam;
                } else {
                    *d2 = *d2 / gam.powi(2);
                    h21 = h21 * gam;
                    h22 = h22 * gam;
                }
            }
        }
    }

    if flag < zero {
        param[[1]] = h11;
        param[[2]] = h21;
        param[[3]] = h12;
        param[[4]] = h22;
    } else if flag == zero {
        param[[2]] = h21;
        param[[3]] = h12;
    } else {
        param[[1]] = h11;
        param[[4]] = h22;
    }

    param[[1]] = flag;
}
