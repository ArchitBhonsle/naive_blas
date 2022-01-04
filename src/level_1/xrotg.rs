use num_traits::Float;

/// Computes the parameters for a Givens rotation.
///
/// Given the co-ordinates (a, b), this function returns the c, s, r and z such that:
///
/// |  c s | | a | = | r |
/// | -s c | | b |   | 0 |
///
/// and z is
///     s   if |a| > |b|,
///     1/c if c != 0,
///     1   otherwise.
///
/// Input:
///     a: x co-ordinate for the point
///     b: y co-ordinate for the point
/// Output:
///     a: contains the parameter r
///     b: contains the parameter z
///     c: contains the parameter c
///     s: contains the parameter s
///
/// Note: Intel MKL also defines rotg for c and z.
pub fn rrotg<T: Float>(a: &mut T, b: &mut T, c: &mut T, s: &mut T) {
    let zero = num_traits::zero();
    let one = num_traits::one();

    let safmin = T::min_value();
    let safmax = T::max_value();

    let anorm = a.abs();
    let bnorm = b.abs();

    if bnorm == zero {
        *c = one;
        *s = zero;
        *b = zero;
    } else if anorm == zero {
        *c = zero;
        *s = one;
        *a = *b;
        *b = one;
    } else {
        let scl = T::min(safmax, T::max(T::max(safmin, anorm), bnorm));

        let sigma = if anorm > bnorm {
            a.signum()
        } else {
            b.signum()
        };

        let r = sigma * (scl * T::sqrt(T::powi(*a / scl, 2) + T::powi(*b / scl, 2)));
        *c = *a / r;
        *s = *b / r;

        let z = if anorm > bnorm {
            *s
        } else if *c != zero {
            one / *c
        } else {
            one
        };

        *a = r;
        *b = z;
    }
}
