/**
Construct the modified Givens transformation matrix `h` which zeros the second
component of the 2-vector `[d1.sqrt(). d2.sqrt()]` with `param[0] = dflag`.

There are four possibilities for `h`:
```
flag    -1.0     0.0       1        -1

h     h11 h12  1.0 h12  d11  1.0  1.0 0.0
      h21 h22  h21 1.0  -1.0 h22  0.0 1.0
```

Note that the constant values are implied from the `param[0]` and not stored in
`param`.

Values for `gamsq` and `rgamsq` used may be inexact since the actual scaling is
done using `gam`.

# Arguments
d1, d2, x1, y1 are mutable references to f32s

# Returns (arguments themselves are used as "return values")
d1, d2, dx1 -> TODO
param -> [flag, h11, h21, h12, h22]

Based on reference BLAS level1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
pub extern "C" fn srotmg(d1: &mut f32, d2: &mut f32, x1: &mut f32, y1: &mut f32, param: *mut f32) {
    let gam = 4096_f32;
    let gamsq = 1.67772e7_f32;
    let rgamsq = 5.96046e-8_f32;

    unsafe {
        let param = std::slice::from_raw_parts_mut(param, 5);

        let mut flag: f32;
        // There might be a bug regarding the following values.
        // Removing intializations here causes "use of possibily-uninitialized" error
        let (mut h11, mut h12, mut h21, mut h22) = (0_f32, 0_f32, 0_f32, 0_f32);

        if *d1 < 0. {
            flag = -1.;
            h11 = 0.;
            h12 = 0.;
            h21 = 0.;
            h22 = 0.;

            *d1 = 0.;
            *d2 = 0.;
            *x1 = 0.;
        } else {
            let p2 = *d2 * *y1;

            if p2 == 0. {
                flag = -2.;
                param[0] = flag;
                return;
            }

            let p1 = *d1 * *x1;
            let q2 = p2 * *y1;
            let q1 = p1 * *x1;

            if q1.abs() > q2.abs() {
                h21 = -*y1 / *x1;
                h12 = p2 / p1;

                let u = 1. - h12 * h21;

                if u > 0. {
                    flag = 0.;
                    *d1 = *d1 / u;
                    *d2 = *d2 / u;
                    *x1 = *x1 * u;
                } else {
                    // Added for safety. Only triggered by to rounding errors
                    // Check DOI: 10.1145/355841.355847
                    flag = -1.;
                    h11 = 0.;
                    h12 = 0.;
                    h21 = 0.;
                    h22 = 0.;

                    *d1 = 0.;
                    *d2 = 0.;
                    *x1 = 0.;
                }
            } else {
                if q2 < 0. {
                    flag = -1.;
                    h11 = 0.;
                    h12 = 0.;
                    h21 = 0.;
                    h22 = 0.;

                    *d1 = 0.;
                    *d2 = 0.;
                    *x1 = 0.;
                } else {
                    flag = 1.;
                    h11 = p1 / p2;
                    h22 = *x1 / *y1;
                    let u = 1. + h11 * h22;
                    let temp = *d1 / u;
                    *d2 = *d1 / u;
                    *d1 = temp;
                    *x1 = *y1 * u;
                }
            }
            // scale check
            if *d1 != 0. {
                while (*d1 < rgamsq) || (*d1 > gamsq) {
                    if flag == 0. {
                        flag = -1.;
                        h11 = 1.;
                        h22 = 1.;
                    } else {
                        flag = -1.;
                        h21 = -1.;
                        h12 = 1.;
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

            if *d2 != 0. {
                while (d2.abs() < rgamsq) || (d2.abs() > gamsq) {
                    if flag == 0. {
                        flag = -1.;
                        h11 = 1.;
                        h22 = 1.;
                    } else {
                        flag = -1.;
                        h21 = -1.;
                        h12 = 1.;
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
        };

        if flag < 0. {
            param[1] = h11;
            param[2] = h21;
            param[3] = h12;
            param[4] = h22;
        } else if flag == 0. {
            param[2] = h12;
            param[3] = h22;
        } else {
            param[1] = h11;
            param[4] = h21;
        }
        param[0] = flag;
    }
}

/**
Construct the modified Givens transformation matrix `h` which zeros the second
component of the 2-vector `[d1.sqrt(). d2.sqrt()]` with `param[0] = dflag`.

There are four possibilities for `h`:
```
flag    -1.0     0.0       1        -1

h     h11 h12  1.0 h12  d11  1.0  1.0 0.0
      h21 h22  h21 1.0  -1.0 h22  0.0 1.0
```

Note that the constant values are implied from the `param[0]` and not stored in
`param`.

Values for `gamsq` and `rgamsq` used may be inexact since the actual scaling is
done using `gam`.

# Arguments
d1, d2, x1, y1 are mutable references to f64s

# Returns (arguments themselves are used as "return values")
d1, d2, dx1 -> TODO
param -> [flag, h11, h21, h12, h22]

Based on reference BLAS level1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
pub extern "C" fn drotmg(d1: &mut f64, d2: &mut f64, x1: &mut f64, y1: &mut f64, param: *mut f64) {
    let gam = 4096_f64;
    let gamsq = 16777216_f64;
    let rgamsq = 5.9604645e-8_f64;

    unsafe {
        let param = std::slice::from_raw_parts_mut(param, 5);

        let mut flag: f64;
        // There might be a bug regarding the following values.
        // Removing intializations here causes "use of possibily-uninitialized" error
        let (mut h11, mut h12, mut h21, mut h22) = (0_f64, 0_f64, 0_f64, 0_f64);

        if *d1 < 0. {
            flag = -1.;
            h11 = 0.;
            h12 = 0.;
            h21 = 0.;
            h22 = 0.;

            *d1 = 0.;
            *d2 = 0.;
            *x1 = 0.;
        } else {
            let p2 = *d2 * *y1;

            if p2 == 0. {
                flag = -2.;
                param[0] = flag;
                return;
            }

            let p1 = *d1 * *x1;
            let q2 = p2 * *y1;
            let q1 = p1 * *x1;

            if q1.abs() > q2.abs() {
                h21 = -*y1 / *x1;
                h12 = p2 / p1;

                let u = 1. - h12 * h21;

                if u > 0. {
                    flag = 0.;
                    *d1 = *d1 / u;
                    *d2 = *d2 / u;
                    *x1 = *x1 * u;
                } else {
                    // Added for safety. Only triggered by to rounding errors
                    // Check DOI: 10.1145/355841.355847
                    flag = -1.;
                    h11 = 0.;
                    h12 = 0.;
                    h21 = 0.;
                    h22 = 0.;

                    *d1 = 0.;
                    *d2 = 0.;
                    *x1 = 0.;
                }
            } else {
                if q2 < 0. {
                    flag = -1.;
                    h11 = 0.;
                    h12 = 0.;
                    h21 = 0.;
                    h22 = 0.;

                    *d1 = 0.;
                    *d2 = 0.;
                    *x1 = 0.;
                } else {
                    flag = 1.;
                    h11 = p1 / p2;
                    h22 = *x1 / *y1;
                    let u = 1. + h11 * h22;
                    let temp = *d1 / u;
                    *d2 = *d1 / u;
                    *d1 = temp;
                    *x1 = *y1 * u;
                }
            }
            // scale check
            if *d1 != 0. {
                while (*d1 < rgamsq) || (*d1 > gamsq) {
                    if flag == 0. {
                        flag = -1.;
                        h11 = 1.;
                        h22 = 1.;
                    } else {
                        flag = -1.;
                        h21 = -1.;
                        h12 = 1.;
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

            if *d2 != 0. {
                while (d2.abs() < rgamsq) || (d2.abs() > gamsq) {
                    if flag == 0. {
                        flag = -1.;
                        h11 = 1.;
                        h22 = 1.;
                    } else {
                        flag = -1.;
                        h21 = -1.;
                        h12 = 1.;
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
        };

        if flag < 0. {
            param[1] = h11;
            param[2] = h21;
            param[3] = h12;
            param[4] = h22;
        } else if flag == 0. {
            param[2] = h12;
            param[3] = h22;
        } else {
            param[1] = h11;
            param[4] = h21;
        }
        param[0] = flag;
    }
}
