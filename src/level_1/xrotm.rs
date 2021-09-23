/**
Applies modified Givens rotation of points in the plane

Given two vectors `x` and `y`, each element of these vectors is replaced as
follows:
```
[ x[i] ] = h [ x[i] ]
[ y[i] ]     [ y[i] ]
```

`n`, `inc_x` and `inc_y` are `isize`s
`x` and `y` are arrays of `f32`s
`param` is an array of 5 `f32`s

# Arguments
`n` -> specifies the number of elements in each array
`x` -> array of size at least `(n-1) * inc_x.abs() + 1`
`inc_x` -> specifies the increments for elements of `x`
`y` -> array of size at least `(n-1) * inc_y.abs() + 1`
`inc_y` -> specifies the increments for elements of `y`
`param` -> describes `h`, `param[0]` is `flag` based on which the elements of
`h` as shown in the table. `h11`, `h21`, `h12`, `h22` are `param[1-4]`. The
constant values are assumed based on the `flag` and may not be reflected in `param`
```
+------+-------------+
| flag |          h  |
+------+-------------+
| -1.0 | [ h11 h12 ] |
|      | [ h21 h22 ] |
|------|-------------|
|  0.0 | [ 1.0 h12 ] |
|      | [ h21 1.0 ] |
|------|-------------|
|  1.0 | [ h11 1.0 ] |
|      | [-1.0 h22 ] |
|------|-------------|
| -2.0 | [ 1.0 0.0 ] |
|      | [ 0.0 1.0 ] |
+------+-------------+
```

# Returns (arguments themselves are used as "return values")
`x` -> `x` where each element is replaced by `h11*x[i] + h12*y[i]`
`y` -> `y` where each element is replaced by `h21*x[i] + h22*y[i]`


Based on reference BLAS level-1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
#[no_mangle]
pub extern "C" fn srotm(
    n: &mut isize,
    x: *mut f32,
    inc_x: &mut isize,
    y: *mut f32,
    inc_y: &mut isize,
    param: *mut f32,
) {
    let (n, inc_x, inc_y) = (*n, *inc_x, *inc_y);

    unsafe {
        let x = std::slice::from_raw_parts_mut(x, n as usize);
        let y = std::slice::from_raw_parts_mut(y, n as usize);
        let param = std::slice::from_raw_parts_mut(param, 5);

        let flag = param[0];
        if n < 0 || flag + 2. == 0. {
            return;
        }

        if inc_x == inc_y && inc_x > 0 {
            // equal and positive increments
            let n_steps = n * inc_x;

            if flag < 0. {
                let (h11, h12, h21, h22) = (param[1], param[3], param[2], param[4]);
                for i in (0..n_steps).step_by(inc_x as usize) {
                    let i = i as usize;
                    let (w, z) = (x[i], y[i]);

                    x[i] = w * h11 + z * h12;
                    y[i] = w * h21 + z * h22;
                }
            } else if flag == 0. {
                let (h12, h21) = (param[3], param[2]);
                for i in (0..n_steps).step_by(inc_x as usize) {
                    let i = i as usize;
                    let (w, z) = (x[i], y[i]);

                    x[i] = w + z * h12;
                    y[i] = w * h21 + z;
                }
            } else {
                let (h11, h22) = (param[1], param[4]);
                for i in (0..n_steps).step_by(inc_x as usize) {
                    let i = i as usize;
                    let (w, z) = (x[i], y[i]);

                    x[i] = w * h11 + z;
                    y[i] = -w + h22 * z;
                }
            }
        } else {
            // all other cases
            let (kx, ky) = (
                if inc_x < 0 { 1 + (1 - n) * inc_x } else { 1 },
                if inc_y < 0 { 1 + (1 - n) * inc_y } else { 1 },
            );
            let (mut kx, mut ky) = (kx as isize, ky as isize);

            if flag < 0. {
                let (h11, h12, h21, h22) = (param[1], param[3], param[2], param[4]);

                for _ in 0..n {
                    let (kxi, kyi) = (kx as usize, ky as usize);
                    let (w, z) = (x[kxi], y[kyi]);

                    x[kxi] = w * h11 + z * h12;
                    y[kyi] = w * h21 + z * h22;

                    kx += inc_x;
                    ky += inc_y
                }
            } else if flag == 0. {
                let (h12, h21) = (param[3], param[2]);

                for _ in 0..n {
                    let (kxi, kyi) = (kx as usize, ky as usize);
                    let (w, z) = (x[kxi], y[kyi]);

                    x[kxi] = w + z * h12;
                    y[kyi] = w * h21 + z;

                    kx += inc_x;
                    ky += inc_y
                }
            } else {
                let (h11, h22) = (param[1], param[4]);

                for _ in 0..n {
                    let (kxi, kyi) = (kx as usize, ky as usize);
                    let (w, z) = (x[kxi], y[kyi]);

                    x[kxi] = w * h11 + z;
                    y[kyi] = -w + z * h22;

                    kx += inc_x;
                    ky += inc_y
                }
            }
        }
    }
}

/**
Applies modified Givens rotation of points in the plane

Given two vectors `x` and `y`, each element of these vectors is replaced as
follows:
```
[ x[i] ] = h [ x[i] ]
[ y[i] ]     [ y[i] ]
```

`n`, `inc_x` and `inc_y` are `isize`s
`x` and `y` are arrays of `f64`s
`param` is an array of 5 `f64`s

# Arguments
`n` -> specifies the number of elements in each array
`x` -> array of size at least `(n-1) * inc_x.abs() + 1`
`inc_x` -> specifies the increments for elements of `x`
`y` -> array of size at least `(n-1) * inc_y.abs() + 1`
`inc_y` -> specifies the increments for elements of `y`
`param` -> describes `h`, `param[0]` is `flag` based on which the elements of
`h` as shown in the table. `h11`, `h21`, `h12`, `h22` are `param[1-4]`. The
constant values are assumed based on the `flag` and may not be reflected in `param`
```
+------+-------------+
| flag |          h  |
+------+-------------+
| -1.0 | [ h11 h12 ] |
|      | [ h21 h22 ] |
|------|-------------|
|  0.0 | [ 1.0 h12 ] |
|      | [ h21 1.0 ] |
|------|-------------|
|  1.0 | [ h11 1.0 ] |
|      | [-1.0 h22 ] |
|------|-------------|
| -2.0 | [ 1.0 0.0 ] |
|      | [ 0.0 1.0 ] |
+------+-------------+
```

# Returns (arguments themselves are used as "return values")
`x` -> `x` where each element is replaced by `h11*x[i] + h12*y[i]`
`y` -> `y` where each element is replaced by `h21*x[i] + h22*y[i]`


Based on reference BLAS level-1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
#[no_mangle]
pub extern "C" fn drotm(
    n: &mut isize,
    x: *mut f64,
    inc_x: &mut isize,
    y: *mut f64,
    inc_y: &mut isize,
    param: *mut f64,
) {
    let (n, inc_x, inc_y) = (*n, *inc_x, *inc_y);

    unsafe {
        let x = std::slice::from_raw_parts_mut(x, n as usize);
        let y = std::slice::from_raw_parts_mut(y, n as usize);
        let param = std::slice::from_raw_parts_mut(param, 5);

        let flag = param[0];
        if n < 0 || flag + 2. == 0. {
            return;
        }

        if inc_x == inc_y && inc_x > 0 {
            // equal and positive increments
            let n_steps = n * inc_x;

            if flag < 0. {
                let (h11, h12, h21, h22) = (param[1], param[3], param[2], param[4]);

                for i in (0..n_steps).step_by(inc_x as usize) {
                    let i = i as usize;
                    let (w, z) = (x[i], y[i]);

                    x[i] = w * h11 + z * h12;
                    y[i] = w * h21 + z * h22;
                }
            } else if flag == 0. {
                let (h12, h21) = (param[3], param[2]);

                for i in (0..n_steps).step_by(inc_x as usize) {
                    let i = i as usize;
                    let (w, z) = (x[i], y[i]);

                    x[i] = w + z * h12;
                    y[i] = w * h21 + z;
                }
            } else {
                let (h11, h22) = (param[1], param[4]);

                for i in (0..n_steps).step_by(inc_x as usize) {
                    let i = i as usize;
                    let (w, z) = (x[i], y[i]);

                    x[i] = w * h11 + z;
                    y[i] = -w + h22 * z;
                }
            }
        } else {
            // all other cases
            let (kx, ky) = (
                if inc_x < 0 { 1 + (1 - n) * inc_x } else { 1 },
                if inc_y < 0 { 1 + (1 - n) * inc_y } else { 1 },
            );
            let (mut kx, mut ky) = (kx as isize, ky as isize);

            if flag < 0. {
                let (h11, h12, h21, h22) = (param[1], param[3], param[2], param[4]);

                for _ in 0..n {
                    let (kxi, kyi) = (kx as usize, ky as usize);
                    let (w, z) = (x[kxi], y[kyi]);

                    x[kxi] = w * h11 + z * h12;
                    y[kyi] = w * h21 + z * h22;

                    kx += inc_x;
                    ky += inc_y
                }
            } else if flag == 0. {
                let (h12, h21) = (param[3], param[2]);

                for _ in 0..n {
                    let (kxi, kyi) = (kx as usize, ky as usize);
                    let (w, z) = (x[kxi], y[kyi]);

                    x[kxi] = w + z * h12;
                    y[kyi] = w * h21 + z;

                    kx += inc_x;
                    ky += inc_y
                }
            } else {
                let (h11, h22) = (param[1], param[4]);

                for _ in 0..n {
                    let (kxi, kyi) = (kx as usize, ky as usize);
                    let (w, z) = (x[kxi], y[kyi]);

                    x[kxi] = w * h11 + z;
                    y[kyi] = -w + z * h22;

                    kx += inc_x;
                    ky += inc_y
                }
            }
        }
    }
}
