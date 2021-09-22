/**
Given the cartesian coordinates `(a, b)` of a point, parameters `c`, `s`, `r`
and `z` are returned such that,
```
[  c  s ] [ a ] = [ r ]
[ -s  c ] [ b ]   [ 0 ]

and

z =   s  , if |a| > |b|
    1/c  , else if c != 0
      z  , else
```

# Argumemts
`a`, `b`, `c`, `s` are mutable references to `f32`s

# Returns (arguments themselves are used as "return values")
`a` -> `r`
`b` -> `z`
`c` -> `c`
`s` -> `s`


Based on reference BLAS level-1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
#[no_mangle]
pub extern "C" fn srotg(a: &mut f32, b: &mut f32, c: &mut f32, s: &mut f32) {
    let a_norm = a.abs();
    let b_norm = b.abs();

    if b_norm == 0. {
        *c = 1.;
        *s = 0.;
        *b = 0.;
        return;
    }

    if a_norm == 0. {
        *c = 0.;
        *s = 1.;
        *a = *b;
        *b = 1.;
        return;
    }

    let scl = f32::MAX.min(f32::MIN.max(a_norm).max(b_norm));
    let sigma = if a_norm > b_norm {
        a.signum()
    } else {
        b.signum()
    };

    let r = sigma * (scl * ((*a / scl).powi(2) + (*b / scl).powi(2)).sqrt());
    *c = *a / r;
    *s = *b / r;

    let z = if a_norm > b_norm {
        *s
    } else if *c != 0. {
        1. / *c
    } else {
        1.
    };

    *a = r;
    *b = z;
}

/**
Given the cartesian coordinates `(a, b)` of a point, parameters `c`, `s`, `r`
and `z` are returned such that,
```
[  c  s ] [ a ] = [ r ]
[ -s  c ] [ b ]   [ 0 ]

and

z =   s  , if |a| > |b|
    1/c  , else if c != 0
      z  , else
```

# Argumemts
`a`, `b`, `c`, `s` are mutable references to `f64`s

# Returns (arguments themselves are used as "return values")
`a` -> `r`
`b` -> `z`
`c` -> `c`
`s` -> `s`


Based on reference BLAS level-1 routine
Reference BLAS is a software package provided by Univ. of Tennessee,
Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd.
**/
#[no_mangle]
pub extern "C" fn drotg(a: &mut f64, b: &mut f64, c: &mut f64, s: &mut f64) {
    let a_norm = a.abs();
    let b_norm = b.abs();

    if b_norm == 0. {
        *c = 1.;
        *s = 0.;
        *b = 0.;
        return;
    }

    if a_norm == 0. {
        *c = 0.;
        *s = 1.;
        *a = *b;
        *b = 1.;
        return;
    }

    let scl = f64::MAX.min(f64::MIN.max(a_norm).max(b_norm));
    let sigma = if a_norm > b_norm {
        a.signum()
    } else {
        b.signum()
    };

    let r = sigma * (scl * ((*a / scl).powi(2) + (*b / scl).powi(2)).sqrt());
    *c = *a / r;
    *s = *b / r;

    let z = if a_norm > b_norm {
        *s
    } else if *c != 0. {
        1. / *c
    } else {
        1.
    };

    *a = r;
    *b = z;
}
