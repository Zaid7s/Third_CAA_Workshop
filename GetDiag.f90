MODULE GetDiag
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: DiagLHS

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         

    SUBROUTINE DiagLHS( iMin        ,& 
                        iMax        ,& 
                        bMin        ,& 
                        bMax        ,& 
                        delta_Zi    ,&
                        delta_t     ,& 
                        LowerDiag   ,& 
                        Diag        ,& 
                        UpperDiag   ,&
                        Scaling_Fac, &
                        A_Jac       ,&
                        S_Jac       ,&
                        Source_Fac  ,&
                        Jac_Curv)

    INTEGER, INTENT(IN)::   iMin         ,&
                            iMax         ,&
                            bMin         ,&
                            bMax
    
    REAL(KIND = rDef), INTENT(IN) ::    Scaling_Fac     ,&
                                        delta_Zi        ,&
                                        delta_t

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Source_Fac       ,&
                                                    Jac_Curv

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    A_Jac      ,&
                                                            S_Jac

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(OUT):: LowerDiag     ,&
                                                          Diag          ,&
                                                          UpperDiag

    INTEGER:: i, j, k

    REAL(KIND = rDef):: epsi

    REAL(KIND = rDef), DIMENSION(iMax, bMax, bMax) ::   Diss_LowerDiag      ,&
                                                        Diss_Diag           ,&
                                                        Diss_UpperDiag      ,&
                                                        Diff_LowerDiag      ,&
                                                        Diff_Diag           ,&
                                                        Diff_UpperDiag

!! The Value of epsi below is dependent on the value of the Dissipation 
!! coefficient from the RHS (check GetDissipation.f90 and GetRHSVector.f90). 
!! epsi_LHS (shown below) must be 2.0*RHS_Coefficnet

    epsi = 0.0_rDef/(delta_t*1024.0_rDef) 

! Main Diagonal
    DO i = iMin, iMax
        IF (i == iMax) THEN
            Diss_Diag(i, 1, 1) = -1.0_rDef*delta_t
            Diss_Diag(i, 1, 2) = -1.0_rDef*delta_t
            Diss_Diag(i, 1, 3) = -1.0_rDef*delta_t
            Diss_Diag(i, 2, 1) = -1.0_rDef*delta_t
            Diss_Diag(i, 2, 2) = -1.0_rDef*delta_t
            Diss_Diag(i, 2, 3) = -1.0_rDef*delta_t
            Diss_Diag(i, 3, 1) = -1.0_rDef*delta_t
            Diss_Diag(i, 3, 2) = -1.0_rDef*delta_t
            Diss_Diag(i, 3, 3) = -1.0_rDef*delta_t
        ELSE
            Diss_Diag(i, 1, 1) = -2.0_rDef*delta_t
            Diss_Diag(i, 1, 2) = -2.0_rDef*delta_t
            Diss_Diag(i, 1, 3) = -2.0_rDef*delta_t
            Diss_Diag(i, 2, 1) = -2.0_rDef*delta_t
            Diss_Diag(i, 2, 2) = -2.0_rDef*delta_t
            Diss_Diag(i, 2, 3) = -2.0_rDef*delta_t
            Diss_Diag(i, 3, 1) = -2.0_rDef*delta_t
            Diss_Diag(i, 3, 2) = -2.0_rDef*delta_t
            Diss_Diag(i, 3, 3) = -2.0_rDef*delta_t
        END IF
    END DO

! Lower Diagonal
    DO i = iMin, iMax
        IF (i == iMin) THEN
            Diss_LowerDiag(i, 1, 1) = 0.0_rDef
            Diss_LowerDiag(i, 1, 2) = 0.0_rDef
            Diss_LowerDiag(i, 1, 3) = 0.0_rDef
            Diss_LowerDiag(i, 2, 1) = 0.0_rDef
            Diss_LowerDiag(i, 2, 2) = 0.0_rDef
            Diss_LowerDiag(i, 2, 3) = 0.0_rDef
            Diss_LowerDiag(i, 3, 1) = 0.0_rDef
            Diss_LowerDiag(i, 3, 2) = 0.0_rDef
            Diss_LowerDiag(i, 3, 3) = 0.0_rDef
        ELSE
            Diss_LowerDiag(i, 1, 1) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 1, 2) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 1, 3) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 2, 1) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 2, 2) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 2, 3) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 3, 1) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 3, 2) = 1.0_rDef*delta_t
            Diss_LowerDiag(i, 3, 3) = 1.0_rDef*delta_t
        END IF
    END DO

! Upper Diagonal
    DO i = iMin, iMax
        IF (i == iMax) THEN
            Diss_UpperDiag(i, 1, 1) = 0.0_rDef
            Diss_UpperDiag(i, 1, 2) = 0.0_rDef
            Diss_UpperDiag(i, 1, 3) = 0.0_rDef
            Diss_UpperDiag(i, 2, 1) = 0.0_rDef
            Diss_UpperDiag(i, 2, 2) = 0.0_rDef
            Diss_UpperDiag(i, 2, 3) = 0.0_rDef
            Diss_UpperDiag(i, 3, 1) = 0.0_rDef
            Diss_UpperDiag(i, 3, 2) = 0.0_rDef
            Diss_UpperDiag(i, 3, 3) = 0.0_rDef
        ELSEIF (i == iMin) THEN
            Diss_UpperDiag(i, 1, 1) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 1, 2) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 1, 3) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 2, 1) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 2, 2) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 2, 3) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 3, 1) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 3, 2) = -1.0_rDef*delta_t
            Diss_UpperDiag(i, 3, 3) = -1.0_rDef*delta_t
        ELSE
            Diss_UpperDiag(i, 1, 1) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 1, 2) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 1, 3) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 2, 1) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 2, 2) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 2, 3) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 3, 1) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 3, 2) = 1.0_rDef*delta_t
            Diss_UpperDiag(i, 3, 3) = 1.0_rDef*delta_t
        END IF
    END DO

! Main Diagonal
    DO i = iMin, iMax
        IF (i == iMin) THEN
            Diff_Diag(i, 1, 1) = 2.0_rDef-delta_t*S_Jac(1, 1, i)*Source_Fac(i)
            Diff_Diag(i, 1, 2) = 0.0_rDef-delta_t*S_Jac(1, 2, i)*Source_Fac(i)
            Diff_Diag(i, 1, 3) = 0.0_rDef-delta_t*S_Jac(1, 3, i)*Source_Fac(i)
            Diff_Diag(i, 2, 1) = 0.0_rDef-delta_t*S_Jac(2, 1, i)*Source_Fac(i)
            Diff_Diag(i, 2, 2) = 2.0_rDef-delta_t*S_Jac(2, 2, i)*Source_Fac(i)
            Diff_Diag(i, 2, 3) = 0.0_rDef-delta_t*S_Jac(2, 3, i)*Source_Fac(i)
            Diff_Diag(i, 3, 1) = 0.0_rDef-delta_t*S_Jac(3, 1, i)*Source_Fac(i)
            Diff_Diag(i, 3, 2) = 0.0_rDef-delta_t*S_Jac(3, 2, i)*Source_Fac(i)
            Diff_Diag(i, 3, 3) = 2.0_rDef-delta_t*S_Jac(3, 3, i)*Source_Fac(i)
        ELSE IF (i == iMax) THEN
            Diff_Diag(i, 1, 1) = (3.0_rDef/2.0_rDef)+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(1, 1, i)&
                                 - delta_t*S_Jac(1, 1, i)*Source_Fac(i)
            Diff_Diag(i, 1, 2) = 0.0_rDef+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(1, 2, i)&
                                 - delta_t*S_Jac(1, 2, i)*Source_Fac(i)
            Diff_Diag(i, 1, 3) = 0.0_rDef+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(1, 3, i)&
                                 - delta_t*S_Jac(1, 3, i)*Source_Fac(i)
            Diff_Diag(i, 2, 1) = 0.0_rDef+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(2, 1, i)&
                                 - delta_t*S_Jac(2, 1, i)*Source_Fac(i)
            Diff_Diag(i, 2, 2) = (3.0_rDef/2.0_rDef)+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(2, 2, i)&
                                 - delta_t*S_Jac(2, 2, i)*Source_Fac(i)
            Diff_Diag(i, 2, 3) = 0.0_rDef+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(2, 3, i)&
                                 - delta_t*S_Jac(2, 3, i)*Source_Fac(i)
            Diff_Diag(i, 3, 1) = 0.0_rDef+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(3, 1, i)&
                                 - delta_t*S_Jac(3, 1, i)*Source_Fac(i)
            Diff_Diag(i, 3, 2) = 0.0_rDef+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(3, 2, i)&
                                 - delta_t*S_Jac(3, 2, i)*Source_Fac(i)
            Diff_Diag(i, 3, 3) = (3.0_rDef/2.0_rDef)+3.0_rDef*Scaling_Fac*(delta_t/  &
                                (2.0_rDef*delta_Zi))*A_Jac(3, 3, i)&
                                 - delta_t*S_Jac(3, 3, i)*Source_Fac(i)
        ELSE
            Diff_Diag(i, 1, 1) = (3.0_rDef/2.0_rDef)-delta_t*S_Jac(1, 1, i)*Source_Fac(i)
            Diff_Diag(i, 1, 2) = 0.0_rDef-delta_t*S_Jac(1, 2, i)*Source_Fac(i)
            Diff_Diag(i, 1, 3) = 0.0_rDef-delta_t*S_Jac(1, 3, i)*Source_Fac(i)
            Diff_Diag(i, 2, 1) = 0.0_rDef-delta_t*S_Jac(2, 1, i)*Source_Fac(i)
            Diff_Diag(i, 2, 2) = (3.0_rDef/2.0_rDef)-delta_t*S_Jac(2, 2, i)*Source_Fac(i)
            Diff_Diag(i, 2, 3) = 0.0_rDef-delta_t*S_Jac(2, 3, i)*Source_Fac(i)
            Diff_Diag(i, 3, 1) = 0.0_rDef-delta_t*S_Jac(3, 1, i)*Source_Fac(i)
            Diff_Diag(i, 3, 2) = 0.0_rDef-delta_t*S_Jac(3, 2, i)*Source_Fac(i)
            Diff_Diag(i, 3, 3) = (3.0_rDef/2.0_rDef)-delta_t*S_Jac(3, 3, i)*Source_Fac(i)
        END IF
    END DO

! Lower Diagonal
    DO i = iMin, iMax
        IF (i == iMin) THEN
            Diff_LowerDiag(i, 1, 1) = 0.0_rDef
            Diff_LowerDiag(i, 1, 2) = 0.0_rDef
            Diff_LowerDiag(i, 1, 3) = 0.0_rDef
            Diff_LowerDiag(i, 2, 1) = 0.0_rDef
            Diff_LowerDiag(i, 2, 2) = 0.0_rDef
            Diff_LowerDiag(i, 2, 3) = 0.0_rDef
            Diff_LowerDiag(i, 3, 1) = 0.0_rDef
            Diff_LowerDiag(i, 3, 2) = 0.0_rDef
            Diff_LowerDiag(i, 3, 3) = 0.0_rDef
        ELSEIF (i == iMax) THEN
            Diff_LowerDiag(i, 1, 1) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(1, 1, i-1)
            Diff_LowerDiag(i, 1, 2) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(1, 2, i-1)
            Diff_LowerDiag(i, 1, 3) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(1, 3, i-1)
            Diff_LowerDiag(i, 2, 1) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(2, 1, i-1)
            Diff_LowerDiag(i, 2, 2) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(2, 2, i-1)
            Diff_LowerDiag(i, 2, 3) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(2, 3, i-1)
            Diff_LowerDiag(i, 3, 1) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(3, 1, i-1)
            Diff_LowerDiag(i, 3, 2) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(3, 2, i-1)
            Diff_LowerDiag(i, 3, 3) = -2.0_rDef*Scaling_Fac*(delta_t/&
                                                delta_Zi)*A_Jac(3, 3, i-1)
        ELSE
            Diff_LowerDiag(i, 1, 1) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(1, 1, i-1)
            Diff_LowerDiag(i, 1, 2) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(1, 2, i-1)
            Diff_LowerDiag(i, 1, 3) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(1, 3, i-1)
            Diff_LowerDiag(i, 2, 1) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(2, 1, i-1)
            Diff_LowerDiag(i, 2, 2) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(2, 2, i-1)
            Diff_LowerDiag(i, 2, 3) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(2, 3, i-1)
            Diff_LowerDiag(i, 3, 1) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(3, 1, i-1)
            Diff_LowerDiag(i, 3, 2) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(3, 2, i-1)
            Diff_LowerDiag(i, 3, 3) = -Scaling_Fac*(delta_t/&
                                        (2.0_rDef*delta_Zi))*A_Jac(3, 3, i-1)
        END IF
    END DO

! Upper Diagonal
    DO i = iMin, iMax
        IF (i == iMax) THEN
            Diff_UpperDiag(i, 1, 1) = 0.0_rDef
            Diff_UpperDiag(i, 1, 2) = 0.0_rDef
            Diff_UpperDiag(i, 1, 3) = 0.0_rDef
            Diff_UpperDiag(i, 2, 1) = 0.0_rDef
            Diff_UpperDiag(i, 2, 2) = 0.0_rDef
            Diff_UpperDiag(i, 2, 3) = 0.0_rDef
            Diff_UpperDiag(i, 3, 1) = 0.0_rDef
            Diff_UpperDiag(i, 3, 2) = 0.0_rDef
            Diff_UpperDiag(i, 3, 3) = 0.0_rDef
        ELSEIF (i == iMin) THEN
            Diff_UpperDiag(i, 1, 1) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(1, 1, i+1)
            Diff_UpperDiag(i, 1, 2) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(1, 2, i+1)
            Diff_UpperDiag(i, 1, 3) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(1, 3, i+1)
            Diff_UpperDiag(i, 2, 1) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(2, 1, i+1)
            Diff_UpperDiag(i, 2, 2) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(2, 2, i+1)
            Diff_UpperDiag(i, 2, 3) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(2, 3, i+1)
            Diff_UpperDiag(i, 3, 1) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(3, 1, i+1)
            Diff_UpperDiag(i, 3, 2) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(3, 2, i+1)
            Diff_UpperDiag(i, 3, 3) = Scaling_Fac*(delta_t/delta_Zi)&
                                            *A_Jac(3, 3, i+1)
        ELSE
            Diff_UpperDiag(i, 1, 1) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(1, 1, i+1)
            Diff_UpperDiag(i, 1, 2) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(1, 2, i+1)
            Diff_UpperDiag(i, 1, 3) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(1, 3, i+1)
            Diff_UpperDiag(i, 2, 1) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(2, 1, i+1)
            Diff_UpperDiag(i, 2, 2) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(2, 2, i+1)
            Diff_UpperDiag(i, 2, 3) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(2, 3, i+1)
            Diff_UpperDiag(i, 3, 1) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(3, 1, i+1)
            Diff_UpperDiag(i, 3, 2) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(3, 2, i+1)
            Diff_UpperDiag(i, 3, 3) = Scaling_Fac*(delta_t/(2.0_rDef*delta_Zi))&
                                        *A_Jac(3, 3, i+1)
        END IF
    END DO

    DO k = bMin, bMax
        DO j = bMin, bMax
            DO i = iMin, iMax
                LowerDiag(i, j, k)  = Diff_LowerDiag(i, j, k) + &
                                      epsi*Diss_LowerDiag(i, j, k)
                Diag(i, j, k)       = Diff_Diag(i, j, k)      + &
                                      epsi*Diss_Diag(i, j, k)
                UpperDiag(i, j, k)  = Diff_UpperDiag(i, j, k) + &
                                      epsi*Diss_UpperDiag(i, j, k)
            END DO
        END DO
    END DO

    RETURN

    END SUBROUTINE DiagLHS

END MODULE GetDiag      
