MODULE GetUMatrix
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: UMatrix

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         

    SUBROUTINE UMatrix( iMin            ,&
                        iMax            ,&
                        bMin            ,&
                        bMax            ,&
                        delta_x         ,&
                        delta_Zi        ,&
                        delta_t         ,&
                        A_Plus          ,&
                        A_Minus         ,&
                        S_Jac           ,&
                        S_Plus          ,&
                        S_Minus         ,&
                        Source_Fac      ,&
                        U_UpperDiag     ,&
                        U_Diag          ,&
                        Jac_Curv        ,&
                        Inv_Jac_Curv    ,&
                        dZidX           ,&
                        dZidt)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin

    REAL(KIND = rDef), INTENT(IN) ::    delta_x     ,&
                                        delta_Zi    ,&
                                        delta_t

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt
    
    REAL(KIND = rDef), INTENT(IN), DIMENSION(:, :, :) ::    A_Plus      ,&
                                                            A_Minus     ,&
                                                            S_Jac       ,&
                                                            S_Plus      ,&
                                                            S_Minus

    REAL(KIND = rDef), INTENT(INOUT), DIMENSION(:, :, :) ::   U_UpperDiag  ,&
                                                            U_Diag

    INTEGER:: i, j

    DO i = iMin, iMax
        U_Diag(1, 1, i) = 1.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(1, 1, i) - A_Minus(1, 1, i)) -     &
                            delta_t*S_Jac(1, 1, i)*Source_Fac(i)*0.0_rDef
        U_Diag(1, 2, i) = 0.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(1, 2, i) - A_Minus(1, 2, i)) -     &
                            delta_t*S_Jac(1, 2, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(1, 3, i) = 0.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(1, 3, i) - A_Minus(1, 3, i)) -     &
                            delta_t*S_Jac(1, 3, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(2, 1, i) = 0.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(2, 1, i) - A_Minus(2, 1, i)) -     &
                            delta_t*S_Jac(2, 1, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(2, 2, i) = 1.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(2, 2, i) - A_Minus(2, 2, i)) -     &
                            delta_t*S_Jac(2, 2, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(2, 3, i) = 0.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(2, 3, i) - A_Minus(2, 3, i)) -     &
                            delta_t*S_Jac(2, 3, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(3, 1, i) = 0.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(3, 1, i) - A_Minus(3, 1, i)) -     &
                            delta_t*S_Jac(3, 1, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(3, 2, i) = 0.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(3, 2, i) - A_Minus(3, 2, i)) -     &
                            delta_t*S_Jac(3, 2, i)*Source_Fac(i)*0.0_rDef   
        U_Diag(3, 3, i) = 1.0_rDef+(delta_t/delta_Zi)*    &
                             (A_Plus(3, 3, i) - A_Minus(3, 3, i)) -     &
                            delta_t*S_Jac(3, 3, i)*Source_Fac(i)*0.0_rDef   
    END DO    
    
    DO i = iMin, iMax
        IF (i == iMax) THEN
            U_UpperDiag(1, 1, i) = 0.0_rDef
            U_UpperDiag(1, 2, i) = 0.0_rDef
            U_UpperDiag(1, 3, i) = 0.0_rDef
            U_UpperDiag(2, 1, i) = 0.0_rDef
            U_UpperDiag(2, 2, i) = 0.0_rDef
            U_UpperDiag(2, 3, i) = 0.0_rDef
            U_UpperDiag(3, 1, i) = 0.0_rDef
            U_UpperDiag(3, 2, i) = 0.0_rDef
            U_UpperDiag(3, 3, i) = 0.0_rDef
        ELSE
            U_UpperDiag(1, 1, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(1, 1, i+1))
            U_UpperDiag(1, 2, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(1, 2, i+1))
            U_UpperDiag(1, 3, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(1, 3, i+1))
            U_UpperDiag(2, 1, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(2, 1, i+1))
            U_UpperDiag(2, 2, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(2, 2, i+1))
            U_UpperDiag(2, 3, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(2, 3, i+1))
            U_UpperDiag(3, 1, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(3, 1, i+1))
            U_UpperDiag(3, 2, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(3, 2, i+1))
            U_UpperDiag(3, 3, i) = &
                                    (delta_t/delta_Zi)*(A_Minus(3, 3, i+1))
        END IF
    END DO

    END SUBROUTINE UMatrix

END MODULE GetUMatrix
