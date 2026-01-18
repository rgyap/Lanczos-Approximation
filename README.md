# Lanczos-Approximation

A command-line Fortran program that calculates the coefficients of Lanczos' approximation of the Gamma function.

## How to use?

Compile the `lanczos_coefficients.f90` file then execute the output on the command line. It will then ask for two numbers: $k$ and then $g$. 
The parameter $k$ can be any positive integer, but since this program only makes use of quadruple-precision numbers, any value of $k$ of 11 or more would result in garbage.
The parameter $g$ can be any nonnegative real number. 

It will then output $k+1$ numbers, say $a'_1, ..., a'_k$ (in the same order of the output), which are computed depending on the values of $k$ and $g$. With that, for any $z\in ℂ$ such that $\Re(z)>-0.5$, we have the following approximation:

$$\Gamma(z+1) \approx \sqrt{2\pi} {\left(z+g+\tfrac{1}{2}\right)}^{z+\tfrac{1}{2}} e^{-\left(z+g+\tfrac{1}{2}\right)}\left(a'_0 + \frac{a'_1}{z+1} + {...} + \frac{a'_k}{z+k} \right)$$

For example, with those numbers, one can create the following implementation of the Gamma function in Fortran, which should be valid for all complex numbers.

```Fortran
! Let Coeffs be the array that contains the outputted values of lanczos_coefficients.f90

RECURSIVE FUNCTION GAMMA_L(z) RESULT(res)
    IMPLICIT NONE
    COMPLEX(qp), INTENT(IN) :: z
    COMPLEX(qp) :: z1, res, x, t, y1, y2
    INTEGER :: i
        
    IF (REAL(z,qp) <= -0.5) THEN
        ! Use the reflection formula to handle complex numbers z
        ! for which Re(z) <= -0.5
        res = PI / (SIN(PI*z) * GAMMA_L(1.0_qp - z)) 
    ELSE 
        z1 = z - (1.0_qp, 0.0_qp)
        x = Coeffs(1) 
        DO i=2, K+1
            x = x + (Coeffs(i) / (z1 + i - 1))
        END DO
        t = z1 + G + 0.5_qp
        y1 = SQ2 * SQPI * (t ** (z1 + 0.5_qp))
        y2 = y1 * EXP(-t)
        res = y2 * x
    END IF
END FUNCTION GAMMA_L
```

For more technical details about how this works, please refer to the PDF file, which should also include a Python implementation.

---

**Note** (17 January 2026). A Java version of the program that also uses arbitrary-precision arithmetic is now included in another branch.
