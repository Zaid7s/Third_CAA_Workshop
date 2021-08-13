MODULE ForwardSweep

    USE,  INTRINSIC :: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: F_Sweep 

    INTEGER, PARAMETER :: rDef = REAL64

    CONTAINS
         
    SUBROUTINE F_Sweep (L_LowerDiag         ,&
                        L_Diag              ,&
                        iMin                ,&
                        iMax                ,&
                        bMin                ,&
                        bMax                ,&
                        RHS                 ,&
                        Delta_Q_star2)

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    L_LowerDiag ,&
                                                            L_Diag

    INTEGER, INTENT(IN) ::  iMin    ,&
                            iMax    ,&
                            bMin    ,&
                            bMax

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) :: RHS

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE :: Inv_L_Diag    
    
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) :: Delta_Q_star2

    INTEGER :: i, j, ii, m, n, l
    REAL(KIND = rDef) :: FacL, DetInv
    
    REAL(KIND = rDef), DIMENSION(:, :,  :), ALLOCATABLE :: LowerFac
    
    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  A,      &
                                                        B,      &
                                                        CC,     &
                                                        DD,     &
                                                        EE

    ALLOCATE(LowerFac(bMax, iMin, iMax)     ,&
             A(bMax, bMax)                  ,&
             B(bMax, bMin)                  ,&
             CC(bMax, bMin)                 ,&
             DD(bMax, bMin)                 ,&
             EE(bMax, bMin)                 ,&
             Inv_L_Diag(bMax, bMax, iMax))

    DO i = iMin, iMax
        LowerFac(1, 1, i) = RHS(i, 1)
        LowerFac(2, 1, i) = RHS(i, 2)
        LowerFac(3, 1, i) = RHS(i, 3)
        DO j = i - 1, iMin, -1
            IF (j == i - 1) THEN
                DO m = iMin, bMax
                    DO n = iMin, bMin
                        DO l = iMin, bMax
                            A(m, l)     = L_LowerDiag(m, l, i)
                            B(l, n)     = RHS(j, l)
                        END DO 
                    END DO
                END DO

                DD = MATMUL(A, B)
                LowerFac(1, 1, i) = LowerFac(1, 1, i) - DD(1, 1)
                LowerFac(2, 1, i) = LowerFac(2, 1, i) - DD(2, 1)  
                LowerFac(3, 1, i) = LowerFac(3, 1, i) - DD(3, 1)  
            ELSE
                LowerFac(1, 1, i) = LowerFac(1, 1, i)
                LowerFac(2, 1, i) = LowerFac(2, 1, i)
                LowerFac(3, 1, i) = LowerFac(3, 1, i)
            END IF
        END DO
       

        DetInv  = 1.0_rDef/(L_Diag(1, 1, i)*L_Diag(2, 2, i)*L_Diag(3, 3, i) - &
                            L_Diag(1, 1, i)*L_Diag(2, 3, i)*L_Diag(3, 2, i) - &
                            L_Diag(1, 2, i)*L_Diag(2, 1, i)*L_Diag(3, 3, i) + &
                            L_Diag(1, 2, i)*L_Diag(2, 3, i)*L_Diag(3, 1, i) + &
                            L_Diag(1, 3, i)*L_Diag(2, 1, i)*L_Diag(3, 2, i) - &
                            L_Diag(1, 3, i)*L_Diag(2, 2, i)*L_Diag(3, 1, i))
        
        Inv_L_Diag(1, 1, i) =  DetInv*(L_Diag(2, 2, i)*L_Diag(3, 3, i) - &
                                       L_Diag(2, 3, i)*L_Diag(3, 2, i)) 
        Inv_L_Diag(2, 1, i) = -DetInv*(L_Diag(2, 1, i)*L_Diag(3, 3, i) - &
                                       L_Diag(2, 3, i)*L_Diag(3, 1, i)) 
        Inv_L_Diag(3, 1, i) =  DetInv*(L_Diag(2, 1, i)*L_Diag(3, 2, i) - &
                                       L_Diag(2, 2, i)*L_Diag(3, 1, i)) 
        Inv_L_Diag(1, 2, i) = -DetInv*(L_Diag(1, 2, i)*L_Diag(3, 3, i) - &
                                       L_Diag(1, 3, i)*L_Diag(3, 2, i)) 
        Inv_L_Diag(2, 2, i) =  DetInv*(L_Diag(1, 1, i)*L_Diag(3, 3, i) - &
                                       L_Diag(1, 3, i)*L_Diag(3, 1, i)) 
        Inv_L_Diag(3, 2, i) = -DetInv*(L_Diag(1, 1, i)*L_Diag(3, 2, i) - &
                                       L_Diag(1, 2, i)*L_Diag(3, 1, i)) 
        Inv_L_Diag(1, 3, i) =  DetInv*(L_Diag(1, 2, i)*L_Diag(2, 3, i) - &
                                       L_Diag(1, 3, i)*L_Diag(2, 2, i)) 
        Inv_L_Diag(2, 3, i) = -DetInv*(L_Diag(1, 1, i)*L_Diag(2, 3, i) - &
                                       L_Diag(1, 3, i)*L_Diag(2, 1, i)) 
        Inv_L_Diag(3, 3, i) =  DetInv*(L_Diag(1, 1, i)*L_Diag(2, 2, i) - &
                                       L_Diag(1, 2, i)*L_Diag(2, 1, i)) 
        
        DO m = iMin, bMax
            DO n = iMin, bMin
                DO l = iMin, bMax
                    A(m, l)     = Inv_L_Diag(m, l, i)
                    B(l, n)     = LowerFac(l, n, i)
                END DO 
            END DO
        END DO

        EE = MATMUL(A, B)

        RHS(i, 1)    = EE(1, 1)
        RHS(i, 2)    = EE(2, 1)
        RHS(i, 3)    = EE(3, 1)
        
        Delta_Q_star2(i, 1) = RHS(i, 1)
        Delta_Q_star2(i, 2) = RHS(i, 2)
        Delta_Q_star2(i, 3) = RHS(i, 3)
    
    END DO

    END SUBROUTINE F_Sweep 

END MODULE ForwardSweep
