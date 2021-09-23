/**
Applies a plane rotation

Given two vectors `x` and `y`, each element of these arrays are replaced as
follows:
```
x[i] = c*x[i] + s*y[i]
y[i] = c*y[i] - s*x[i]
```

`n`, `inc_x` and `inc_y` are `i32`
`c` and `s` are `f32`s
`x` and `y` are arrays of `f32`s

# Arguments
`n` -> number of elements in `x` and `y`
`x` -> array, size at least `(n-1) * inc_x.abs() + 1`
`inc_x` -> specifies the increment for `x`
`y` -> array, size at least `(n-1) * inc_y.abs() + 1`
`inc_y` -> specifies the increment for `y`
`c` -> one component defining the rotation
`s` -> one component defining the rotation

# Returns (arguments themselves are used as "return values")
`x` -> `x` after performing the rotation
`y` -> `y` after performing the rotation


Based on reference BLAS level-1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
#[no_mangle]
pub extern "C" fn srotg(
    n: &mut isize,
    x: *mut f32,
    inc_x: &mut isize,
    y: *mut f32,
    inc_y: &mut isize,
    c: &mut f32,
    s: &mut f32,
) {
    let (n, inc_x, inc_y, c, s) = (*n, *inc_x, *inc_y, *c, *s);

    unsafe {
        let x = std::slice::from_raw_parts_mut(x, n as usize);
        let y = std::slice::from_raw_parts_mut(y, n as usize);

        if n < 0 {
            return;
        }

        // When both increments are 1
        if inc_x == 1 && inc_y == 1 {
            for i in 0..(n as usize) {
                let temp = c * x[i] + s * y[i];
                y[i] = c * y[i] - s * x[i];
                x[i] = temp;
            }

            return;
        }

        // All other cases
        let (mut i_x, mut i_y) = (
            if inc_x < 0 { (1 - n) * inc_x + 1 } else { 1 },
            if inc_y < 0 { (1 - n) * inc_y + 1 } else { 1 },
        );

        for _ in 0..(n as usize) {
            let (i_xu, i_yu) = (i_x as usize, i_y as usize);

            let temp = c * x[i_xu] + s * y[i_yu];
            y[i_yu] = c * y[i_yu] - s * x[i_xu];
            x[i_xu] = temp;

            i_x += inc_x;
            i_y += inc_y;
        }
    }
}

/**
Applies a plane rotation

Given two vectors `x` and `y`, each element of these arrays are replaced as
follows:
```
x[i] = c*x[i] + s*y[i]
y[i] = c*y[i] - s*x[i]
```

`n`, `inc_x` and `inc_y` are `isize`
`c` and `s` are `f64`s
`x` and `y` are arrays of `f64`s

# Arguments
`n` -> number of elements in `x` and `y`
`x` -> array, size at least `(n-1) * inc_x.abs() + 1`
`inc_x` -> specifies the increment for `x`
`y` -> array, size at least `(n-1) * inc_y.abs() + 1`
`inc_y` -> specifies the increment for `y`
`c` -> one component defining the rotation
`s` -> one component defining the rotation

# Returns (arguments themselves are used as "return values")
`x` -> `x` after performing the rotation
`y` -> `y` after performing the rotation


Based on reference BLAS level-1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
#[no_mangle]
pub extern "C" fn drotg(
    n: &mut isize,
    x: *mut f64,
    inc_x: &mut isize,
    y: *mut f64,
    inc_y: &mut isize,
    c: &mut f64,
    s: &mut f64,
) {
    let (n, inc_x, inc_y, c, s) = (*n, *inc_x, *inc_y, *c, *s);

    unsafe {
        let x = std::slice::from_raw_parts_mut(x, n as usize);
        let y = std::slice::from_raw_parts_mut(y, n as usize);

        if n < 0 {
            return;
        }

        // When both increments are 1
        if inc_x == 1 && inc_y == 1 {
            for i in 0..(n as usize) {
                let temp = c * x[i] + s * y[i];
                y[i] = c * y[i] - s * x[i];
                x[i] = temp;
            }

            return;
        }

        // All other cases
        let (mut i_x, mut i_y) = (
            if inc_x < 0 { (1 - n) * inc_x + 1 } else { 1 },
            if inc_y < 0 { (1 - n) * inc_y + 1 } else { 1 },
        );

        for _ in 0..(n as usize) {
            let (i_xu, i_yu) = (i_x as usize, i_y as usize);

            let temp = c * x[i_xu] + s * y[i_yu];
            y[i_yu] = c * y[i_yu] - s * x[i_xu];
            x[i_xu] = temp;

            i_x += inc_x;
            i_y += inc_y;
        }
    }
}
