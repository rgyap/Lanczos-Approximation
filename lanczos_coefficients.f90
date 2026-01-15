!  lanczos_coefficients.f90 
!

!****************************************************************************
!
!  PROGRAM: lanczos_coefficients
!
!  PURPOSE:  Calculate the coefficients of Lanczos' approximation of the Gamma function.
!
!****************************************************************************

PROGRAM lanczos_coefficients
     
    IMPLICIT NONE
    INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(p=33, r=4931)
        
    ! Variables
        
    INTEGER :: k
    REAL(qp) :: g
    
    ! Constants
    
    REAL(qp), PARAMETER :: SQ2 = 1.4142135623730950488016887242096980785696718753769480731766797379_qp
    REAL(qp), PARAMETER :: PI = 3.1415926535897932384626433832795028841971693993751058209749445923_qp
    REAL(qp), PARAMETER :: SQPI = 1.7724538509055160272981674833411551827975494561223871282138077898_qp
   
    ! User Input
    WRITE(*,'(A)',advance='no') "Enter value of positive integer parameter 'k': "
    READ(*,*) k
    WRITE(*,'(A)',advance='no') "Enter value of nonnegative real parameter 'g': "
    READ(*,*) g
 
    ! Output
    WRITE(*,*) LANCZOS_COEFFS(k, g)
    
    CONTAINS            
    
    ! Calculations to get the vector "f"
    
    FUNCTION GAMMA_PLUS_HALF(n) RESULT(res)  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        INTEGER :: temp
        REAL(qp) :: res, temp2, n2
        INTEGER :: q
        
        temp = 1
        DO q = 1, n
            temp = temp * (2*q - 1)
        END DO
        
        temp2 = REAL(temp,qp)
        n2 = REAL(n,qp)
        
        res = temp2 * (SQPI / (2.0_qp ** N2))
    END FUNCTION GAMMA_PLUS_HALF
    
    FUNCTION HELPER(a, g) RESULT(res)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: a
        REAL(qp), INTENT(IN) :: g
        REAL(qp) :: a2, res, t0, t1, t2, t3
        
        a2 = REAL(a,qp)
        t0 = SQ2 / PI
        t1 = (a2 + g + 0.5_qp)**(-a2 - 0.5_qp)
        t2 = EXP(a2 + g + 0.5_qp)
        t3 = GAMMA_PLUS_HALF(a)

        res = t0 * t1 * t2 * t3
    END FUNCTION HELPER
    
    FUNCTION VECTOR_f(k,g) RESULT(f)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: k
        REAL(qp), INTENT(IN) :: g
        REAL(qp), DIMENSION(k+1,1) :: f
        INTEGER :: i
        
        DO i = 0, k
            f(i+1,1) = HELPER(i,g)
        END DO
    END FUNCTION VECTOR_f
    
    ! Calculations to get the matrix "C"
    
    FUNCTION CHEBYSHEV_COEFFS(k) RESULT(dp)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: k
        REAL(qp), DIMENSION(2*k+1,2*k+1) :: dp
        INTEGER :: i, j
        
        dp = 0
        dp(1,1) = 1.0_qp
        dp(2,2) = 1.0_qp
         
        DO i=3, 2*k+1
            dp(i,i) = 2.0_qp * dp(i-1,i-1)
            dp(i,1) = -1.0_qp * dp(i-2,1)
        END DO
       
        DO j=2, 2*k+1
            DO i=j+1, 2*k+1
                dp(i,j) = 2.0_qp * dp(i-1,j-1) - dp(i-2,j)
            END DO 
        END DO 
    END FUNCTION CHEBYSHEV_COEFFS
    
    FUNCTION MATRIX_C(k) RESULT(C)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: k
        REAL(qp), DIMENSION(2*k+1,2*k+1) :: MATRIX
        REAL(qp), DIMENSION(k+1,k+1) :: C
        INTEGER :: i,j
        
        MATRIX = CHEBYSHEV_COEFFS(k)
        
        DO i=1, k+1
            DO j=1, k+1
                C(i,j) = MATRIX(2*i-1,2*j-1)
            END DO
        END DO
        
        C(1,1) = 0.5_qp
        
    END FUNCTION MATRIX_C
    
    ! Calculations to get the Matrix B
    
    RECURSIVE FUNCTION FACTORIAL(n) RESULT(res)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL(qp) :: res
        
        IF (n <= 1) THEN
            res = 1.0_qp
        ELSE
            res = REAL(n,qp) * FACTORIAL(n-1)
        END IF
    END FUNCTION FACTORIAL
    
    FUNCTION MATRIX_B(k) RESULT(B)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: k
        REAL(qp), DIMENSION(k+1,k+1) :: B
        REAL(qp) :: v1, v2, d1
        INTEGER :: x, y, d
        
        B = 0
        
        OUTER: do x=0, k
            INNER: do y=0, k
                IF (x > y) THEN
                    B(x+1,y+1) = 0.0_qp
                    CYCLE INNER
                END IF
                
                IF (x == 0) THEN
                    B(x+1,y+1) = 1.0_qp
                    CYCLE INNER
                END IF
                
                IF (y >= x) THEN
                    d = y - x
                    d1 = REAL(d,qp)
                    v1 = FACTORIAL(y + x - 1) * ((-1.0_qp)**d1)
                    v2 = FACTORIAL(2 * x - 1) * FACTORIAL(d)
                    B(x+1,y+1) = v1 / v2
                    CYCLE INNER
                END IF
            END DO INNER
        END DO OUTER
    END FUNCTION MATRIX_B
    
    ! Calculations to get the Matrix D
    
    FUNCTION MATRIX_D(k) RESULT(D)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: k
        REAL(qp), DIMENSION(k+1,k+1) :: D
        INTEGER :: x, y, x1
        REAL(qp) :: v1, v2
        
        D = 0
        
        OUTER: DO x=0, k
            INNER: DO y=0, k
                IF (x /= y) THEN
                    D(x+1,y+1) = 0.0_qp
                    CYCLE INNER
                END IF
                
                IF (x == 0) THEN
                    D(x+1,y+1) = 1.0_qp
                    CYCLE INNER
                END IF
                
                x1 = x-1
                v1 = -1.0_qp * FACTORIAL(2 * x1 + 1)
                v2 = FACTORIAL(x1) ** 2
                D(x+1,y+1) = v1 / v2
            END DO INNER
        END DO OUTER
    END FUNCTION MATRIX_D
    
    FUNCTION LANCZOS_COEFFS(k, g) RESULT(Coeffs)
        INTEGER, INTENT(IN) :: k
        REAL(qp), INTENT(IN) :: g
        REAL(qp), DIMENSION(k+1,k+1) :: M1, M2
        REAL(qp), DIMENSION(k+1,1) :: f, Coeffs
 
        f = VECTOR_f(k,g)
        M1 = MATMUL(MATRIX_D(k),MATRIX_B(k))
        M2 = MATMUL(M1, MATRIX_C(k))
        Coeffs = MATMUL(M2, f)
    END FUNCTION LANCZOS_COEFFS
            
END PROGRAM lanczos_coefficients

