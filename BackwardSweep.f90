MODULE BackwardSweep

    USE,  INTRINSIC :: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: B_Sweep 

    INTEGER, PARAMETER :: rDef = REAL64

    CONTAINS
         
    SUBROUTINE B_Sweep (U_UpperDiag         ,&
                        U_Diag              ,&
                        iMin                ,&
                        iMax                ,&
                        bMin                ,&
                        bMax                ,&
                        Delta_Q_star        ,&
                        Delta_Q)

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    U_UpperDiag ,&
                                                            U_Diag

    INTEGER, INTENT(IN) ::  iMin    ,&
                            iMax    ,&
                            bMin    ,&
                            bMax

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) :: Delta_Q_star

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE :: Inv_U_Diag    
    
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) :: Delta_Q

    INTEGER :: i, j, ii, m, n, l
    REAL(KIND = rDef) :: FacL, DetInv
    
    REAL(KIND = rDef), DIMENSION(:, :,  :), ALLOCATABLE :: UpperFac
    
    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  A,      &
                                                        B,      &
                                                        CC,     &
                                                        DD,     &
                                                        EE

    ALLOCATE(UpperFac(bMax, iMin, iMax)     ,&
             A(bMax, bMax)                  ,&
             B(bMax, bMin)                  ,&
             CC(bMax, bMin)                 ,&
             DD(bMax, bMin)                 ,&
             EE(bMax, bMin)                 ,&
             Inv_U_Diag(bMax, bMax, iMax))

    DO i = iMax, iMin, -iMin
        UpperFac(1, 1, i) = Delta_Q_star(i, 1)
        UpperFac(2, 1, i) = Delta_Q_star(i, 2)
        UpperFac(3, 1, i) = Delta_Q_star(i, 3)
        DO j = i + 1, iMax
            IF (j == i + 1) THEN
                DO m = iMin, bMax
                    DO n = iMin, bMin
                        DO l = iMin, bMax
                            A(m, l)     = U_UpperDiag(m, l, i)
                            B(l, n)     = Delta_Q_star(j, l)
                        END DO 
                    END DO
                END DO
                DD = MATMUL(A, B)
                UpperFac(1, 1, i) = UpperFac(1, 1, i) - DD(1, 1)
                UpperFac(2, 1, i) = UpperFac(2, 1, i) - DD(2, 1)  
                UpperFac(3, 1, i) = UpperFac(3, 1, i) - DD(3, 1)  
            ELSE
                UpperFac(1, 1, i) = UpperFac(1, 1, i)
                UpperFac(2, 1, i) = UpperFac(2, 1, i)
                UpperFac(3, 1, i) = UpperFac(3, 1, i)
            END IF
        END DO
       
        DetInv  = 1.0_rDef/(U_Diag(1, 1, i)*U_Diag(2, 2, i)*U_Diag(3, 3, i) - &
                            U_Diag(1, 1, i)*U_Diag(2, 3, i)*U_Diag(3, 2, i) - &
                            U_Diag(1, 2, i)*U_Diag(2, 1, i)*U_Diag(3, 3, i) + &
                            U_Diag(1, 2, i)*U_Diag(2, 3, i)*U_Diag(3, 1, i) + &
                            U_Diag(1, 3, i)*U_Diag(2, 1, i)*U_Diag(3, 2, i) - &
                            U_Diag(1, 3, i)*U_Diag(2, 2, i)*U_Diag(3, 1, i))
        
        Inv_U_Diag(1, 1, i) =  DetInv*(U_Diag(2, 2, i)*U_Diag(3, 3, i) - &
                                       U_Diag(2, 3, i)*U_Diag(3, 2, i)) 
        Inv_U_Diag(2, 1, i) = -DetInv*(U_Diag(2, 1, i)*U_Diag(3, 3, i) - &
                                       U_Diag(2, 3, i)*U_Diag(3, 1, i)) 
        Inv_U_Diag(3, 1, i) =  DetInv*(U_Diag(2, 1, i)*U_Diag(3, 2, i) - &
                                       U_Diag(2, 2, i)*U_Diag(3, 1, i)) 
        Inv_U_Diag(1, 2, i) = -DetInv*(U_Diag(1, 2, i)*U_Diag(3, 3, i) - &
                                       U_Diag(1, 3, i)*U_Diag(3, 2, i)) 
        Inv_U_Diag(2, 2, i) =  DetInv*(U_Diag(1, 1, i)*U_Diag(3, 3, i) - &
                                       U_Diag(1, 3, i)*U_Diag(3, 1, i)) 
        Inv_U_Diag(3, 2, i) = -DetInv*(U_Diag(1, 1, i)*U_Diag(3, 2, i) - &
                                       U_Diag(1, 2, i)*U_Diag(3, 1, i)) 
        Inv_U_Diag(1, 3, i) =  DetInv*(U_Diag(1, 2, i)*U_Diag(2, 3, i) - &
                                       U_Diag(1, 3, i)*U_Diag(2, 2, i)) 
        Inv_U_Diag(2, 3, i) = -DetInv*(U_Diag(1, 1, i)*U_Diag(2, 3, i) - &
                                       U_Diag(1, 3, i)*U_Diag(2, 1, i)) 
        Inv_U_Diag(3, 3, i) =  DetInv*(U_Diag(1, 1, i)*U_Diag(2, 2, i) - &
                                       U_Diag(1, 2, i)*U_Diag(2, 1, i)) 
        DO m = iMin, bMax
            DO n = iMin, bMin
                DO l = iMin, bMax
                    A(m, l)     = Inv_U_Diag(m, l, i)
                    B(l, n)     = UpperFac(l, n, i)
                END DO 
            END DO
        END DO

        EE = MATMUL(A, B)

        Delta_Q_star(i, 1)    = EE(1, 1)
        Delta_Q_star(i, 2)    = EE(2, 1)
        Delta_Q_star(i, 3)    = EE(3, 1)
        
        Delta_Q(i, 1) = Delta_Q_star(i, 1)
        Delta_Q(i, 2) = Delta_Q_star(i, 2)
        Delta_Q(i, 3) = Delta_Q_star(i, 3)
    
    END DO

    END SUBROUTINE B_Sweep 

END MODULE BackwardSweep
