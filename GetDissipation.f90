MODULE GetDissipation
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE

    PUBLIC:: Dissipation_RHS

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         
    SUBROUTINE Dissipation_RHS( iMin        ,&
                                iMax        ,&
                                bMin        ,&
                                bMax        ,&
                                delta_X     ,&
                                delta_tau   ,&
                                nD          ,&
                                Dis         ,&
                                Q_n)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           nD

    REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                     delta_tau
                                        
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: Q_n
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: Dis
    
    INTEGER:: i, j
    REAL(KIND = rDef)  :: KCFac
    
    DO j = bMin, bMax
        DO i = iMin, iMax
            Dis(i, j) = 0.0_rDef 
        END DO
    END DO

    IF (nD == 1) THEN  ! Second Order Dissipation, D2
        KCFac = 1.0_rDef/(4.0_rDef)
        DO j = bMin, bMax
            DO i = iMin+1, iMax-1
                Dis(i, j) = KCFac*(1.0_rDef*Q_n(i-1, j)     & 
                                -  2.0_rDef*Q_n(i, j)         &
                                +  1.0_rDef*Q_n(i+1, j))
            END DO
        END DO
        DO j = bMin, bMax
            i = iMin
            Dis(i, j) = KCFac*(  -1.0_rDef*Q_n(i, j)        &
                                + 1.0_rDef*Q_n(i+1, j))
            i = iMax
            Dis(i, j) = KCFac*(  -1.0_rDef*Q_n(i, j)        &
                                + 1.0_rDef*Q_n(i-1, j))
        END DO
    ELSE IF (nD == 2) THEN  !Fourth Order Dissipation, D4
        KCFac = 1.0_rDef/(16.0_rDef)
        DO j = bMin, bMax
            DO i = iMin+2, iMax-2
                Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i-2, j)    & 
                                 +  4.0_rDef*Q_n(i-1, j)    &
                                 -  6.0_rDef*Q_n(i    , j)    &
                                 +  4.0_rDef*Q_n(i+1, j)    &
                                 -  1.0_rDef*Q_n(i+2, j))
            END DO
        END DO
        DO j = bMin, bMax
            
            i = iMin
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  2.0_rDef*Q_n(i+1, j)    &
                             -  1.0_rDef*Q_n(i+2, j))
            
            i = iMin+1
            Dis(i, j) = KCFac*( 2.0_rDef*Q_n(i-1, j)    & 
                             -  5.0_rDef*Q_n(i    , j)    &
                             +  4.0_rDef*Q_n(i+1, j)    &
                             -  1.0_rDef*Q_n(i+2, j))
            
            i = iMax
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  2.0_rDef*Q_n(i-1, j)    &
                             -  1.0_rDef*Q_n(i-2, j))
            
            i = iMax-1
            Dis(i, j) = KCFac*( 2.0_rDef*Q_n(i+1, j)    & 
                             -  5.0_rDef*Q_n(i    , j)    &
                             +  4.0_rDef*Q_n(i-1, j)    &
                             -  1.0_rDef*Q_n(i-2, j))
       END DO
    
    ELSE IF (nD == 3) THEN  !Sixth Order Dissipation, D6
        KCFac = 1.0_rDef/(64.0_rDef)
        DO j = bMin, bMax
            DO i = iMin+3, iMax-3
                Dis(i, j) = KCFac*( 1.0_rDef*Q_n(i-3, j)    & 
                                 -  6.0_rDef*Q_n(i-2, j)    &
                                 + 15.0_rDef*Q_n(i-1, j)    &
                                 - 20.0_rDef*Q_n(i    , j)    &
                                 + 15.0_rDef*Q_n(i+1, j)    &
                                 -  6.0_rDef*Q_n(i+2, j)    &
                                 +  1.0_rDef*Q_n(i+3, j))
            END DO
        END DO
        DO j = bMin, bMax
            
            i = iMin
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  3.0_rDef*Q_n(i+1, j)    &
                             -  3.0_rDef*Q_n(i+2, j)    &
                             +  1.0_rDef*Q_n(i+3, j))
            
            i = iMin+1
            Dis(i, j) = KCFac*( 3.0_rDef*Q_n(i-1, j)    & 
                             - 10.0_rDef*Q_n(i    , j)    &
                             + 12.0_rDef*Q_n(i+1, j)    &
                             -  6.0_rDef*Q_n(i+2, j)    &
                             +  1.0_rDef*Q_n(i+3, j))
            
            i = iMin+2
            Dis(i, j) = KCFac*(-3.0_rDef*Q_n(i-2, j)    & 
                             + 12.0_rDef*Q_n(i-1, j)    &
                             - 19.0_rDef*Q_n(i    , j)    &
                             + 15.0_rDef*Q_n(i+1, j)    &
                             -  6.0_rDef*Q_n(i+2, j)    &
                             +  1.0_rDef*Q_n(i+3, j))
            
            i = iMax
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  3.0_rDef*Q_n(i-1, j)    &
                             -  3.0_rDef*Q_n(i-2, j)    &
                             +  1.0_rDef*Q_n(i-3, j))
            
            i = iMax-1
            Dis(i, j) = KCFac*( 3.0_rDef*Q_n(i+1, j)    & 
                             - 10.0_rDef*Q_n(i    , j)    &
                             + 12.0_rDef*Q_n(i-1, j)    &
                             -  6.0_rDef*Q_n(i-2, j)    &
                             +  1.0_rDef*Q_n(i-3, j))

            i = iMax-2
            Dis(i, j) = KCFac*(-3.0_rDef*Q_n(i+2, j)    & 
                             + 12.0_rDef*Q_n(i+1, j)    &
                             - 19.0_rDef*Q_n(i    , j)    &
                             + 15.0_rDef*Q_n(i-1, j)    &
                             -  6.0_rDef*Q_n(i-2, j)    &
                             +  1.0_rDef*Q_n(i-3, j))
       END DO
    ELSE IF (nD == 4) THEN  !Eighteth Order Dissipation, D8
        KCFac = 1.0_rDef/(256.0_rDef)
        DO j = bMin, bMax
            DO i = iMin+4, iMax-4
                Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i-4, j)    & 
                                 +  8.0_rDef*Q_n(i-3, j)    &
                                 - 28.0_rDef*Q_n(i-2, j)    &
                                 + 56.0_rDef*Q_n(i-1, j)    &
                                 - 70.0_rDef*Q_n(i    , j)    &
                                 + 56.0_rDef*Q_n(i+1, j)    &
                                 - 28.0_rDef*Q_n(i+2, j)    &
                                 +  8.0_rDef*Q_n(i+3, j)    &
                                 -  1.0_rDef*Q_n(i+4, j))
            END DO
        END DO
        DO j = bMin, bMax
            i = iMin
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  4.0_rDef*Q_n(i+1, j)    &
                             -  6.0_rDef*Q_n(i+2, j)    &
                             +  4.0_rDef*Q_n(i+3, j)    &
                             -  1.0_rDef*Q_n(i+4, j))

            i = iMin+1
            Dis(i, j) = KCFac*( 4.0_rDef*Q_n(i-1, j)    & 
                             - 17.0_rDef*Q_n(i    , j)    &
                             + 28.0_rDef*Q_n(i+1, j)    &
                             - 22.0_rDef*Q_n(i+2, j)    &
                             +  8.0_rDef*Q_n(i+3, j)    &
                             -  1.0_rDef*Q_n(i+3, j))

            i = iMin+2
            Dis(i, j) = KCFac*(-6.0_rDef*Q_n(i-2, j)    & 
                             + 28.0_rDef*Q_n(i-1, j)    &
                             - 53.0_rDef*Q_n(i    , j)    &
                             + 52.0_rDef*Q_n(i+1, j)    &
                             - 28.0_rDef*Q_n(i+2, j)    &
                             +  8.0_rDef*Q_n(i+3, j)    &
                             -  1.0_rDef*Q_n(i+4, j))

            i = iMin+3
            Dis(i, j) = KCFac*( 4.0_rDef*Q_n(i-3, j)    & 
                             - 22.0_rDef*Q_n(i-2, j)    &
                             + 52.0_rDef*Q_n(i-1, j)    &
                             - 69.0_rDef*Q_n(i    , j)    &
                             + 56.0_rDef*Q_n(i+1, j)    &
                             - 28.0_rDef*Q_n(i+2, j)    &
                             +  8.0_rDef*Q_n(i+3, j)    &
                             -  1.0_rDef*Q_n(i+4, j))

            i = iMax
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  4.0_rDef*Q_n(i-1, j)    &
                             -  6.0_rDef*Q_n(i-2, j)    &
                             +  4.0_rDef*Q_n(i-3, j)    &
                             -  1.0_rDef*Q_n(i-4, j))

            i = iMax-1
            Dis(i, j) = KCFac*( 4.0_rDef*Q_n(i+1, j)    & 
                             - 17.0_rDef*Q_n(i    , j)    &
                             + 28.0_rDef*Q_n(i-1, j)    &
                             - 22.0_rDef*Q_n(i-2, j)    &
                             +  8.0_rDef*Q_n(i-3, j)    &
                             -  1.0_rDef*Q_n(i-4, j))

            i = iMax-2
            Dis(i, j) = KCFac*(-6.0_rDef*Q_n(i+2, j)    & 
                             + 28.0_rDef*Q_n(i+1, j)    &
                             - 53.0_rDef*Q_n(i    , j)    &
                             + 52.0_rDef*Q_n(i-1, j)    &
                             - 28.0_rDef*Q_n(i-2, j)    &
                             +  8.0_rDef*Q_n(i-3, j)    &
                             -  1.0_rDef*Q_n(i-4, j))

            i = iMax-3
            Dis(i, j) = KCFac*( 4.0_rDef*Q_n(i+3, j)    & 
                             - 22.0_rDef*Q_n(i+2, j)    &
                             + 52.0_rDef*Q_n(i+1, j)    &
                             - 69.0_rDef*Q_n(i    , j)    &
                             + 56.0_rDef*Q_n(i-1, j)    &
                             - 28.0_rDef*Q_n(i-2, j)    &
                             +  8.0_rDef*Q_n(i-3, j)    &
                             -  1.0_rDef*Q_n(i-4, j))
       END DO
    ELSE IF (nD == 5) THEN  !Tenth Order Dissipation, D10
        KCFac = 1.0_rDef/(1024.0_rDef)
        DO j = bMin, bMax
            DO i = iMin+5, iMax-5
                Dis(i, j) = KCFac*(  1.0_rDef*Q_n(i-5, j)    & 
                                 -  10.0_rDef*Q_n(i-4, j)    &
                                 +  45.0_rDef*Q_n(i-3, j)    &
                                 - 120.0_rDef*Q_n(i-2, j)    &
                                 + 210.0_rDef*Q_n(i-1, j)    &
                                 - 252.0_rDef*Q_n(i    , j)    &
                                 + 210.0_rDef*Q_n(i+1, j)    &
                                 - 120.0_rDef*Q_n(i+2, j)    &
                                 +  45.0_rDef*Q_n(i+3, j)    &
                                 -  10.0_rDef*Q_n(i+4, j)    &
                                 +   1.0_rDef*Q_n(i+5, j))
            END DO
        END DO
        DO j = bMin, bMax
            i = iMin
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  5.0_rDef*Q_n(i+1, j)    &
                             - 10.0_rDef*Q_n(i+2, j)    &
                             + 10.0_rDef*Q_n(i+3, j)    &
                             -  5.0_rDef*Q_n(i+4, j)    &
                             +  1.0_rDef*Q_n(i+5, j))

            i = iMin+1
            Dis(i, j) = KCFac*( 5.0_rDef*Q_n(i-1, j)    & 
                             - 26.0_rDef*Q_n(i    , j)    &
                             + 55.0_rDef*Q_n(i+1, j)    &
                             - 60.0_rDef*Q_n(i+2, j)    &
                             + 35.0_rDef*Q_n(i+3, j)    &
                             - 10.0_rDef*Q_n(i+4, j)    &
                             +  1.0_rDef*Q_n(i+5, j))

            i = iMin+2
            Dis(i, j) = KCFac*(-10.0_rDef*Q_n(i-2, j)    & 
                             +  55.0_rDef*Q_n(i-1, j)    &
                             - 126.0_rDef*Q_n(i    , j)    &
                             + 155.0_rDef*Q_n(i+1, j)    &
                             - 110.0_rDef*Q_n(i+2, j)    &
                             +  45.0_rDef*Q_n(i+3, j)    &
                             -  10.0_rDef*Q_n(i+4, j)    &
                             +   1.0_rDef*Q_n(i+5, j))

            i = iMin+3
            Dis(i, j) = KCFac*( 10.0_rDef*Q_n(i-3, j)    & 
                             -  60.0_rDef*Q_n(i-2, j)    &
                             + 155.0_rDef*Q_n(i-1, j)    &
                             - 226.0_rDef*Q_n(i    , j)    &
                             + 205.0_rDef*Q_n(i+1, j)    &
                             - 120.0_rDef*Q_n(i+2, j)    &
                             +  45.0_rDef*Q_n(i+3, j)    &
                             -  10.0_rDef*Q_n(i+4, j)    &
                             +   1.0_rDef*Q_n(i+5, j))

            i = iMin+4
            Dis(i, j) = KCFac*( -5.0_rDef*Q_n(i-4, j)    & 
                             +  35.0_rDef*Q_n(i-3, j)    &
                             - 110.0_rDef*Q_n(i-2, j)    &
                             + 205.0_rDef*Q_n(i-1, j)    &
                             - 251.0_rDef*Q_n(i    , j)    &
                             + 210.0_rDef*Q_n(i+1, j)    &
                             - 120.0_rDef*Q_n(i+2, j)    &
                             +  45.0_rDef*Q_n(i+3, j)    &
                             -  10.0_rDef*Q_n(i+4, j)    &
                             +   1.0_rDef*Q_n(i+5, j))
            
            i = iMax
            Dis(i, j) = KCFac*(-1.0_rDef*Q_n(i    , j)    & 
                             +  5.0_rDef*Q_n(i-1, j)    &
                             - 10.0_rDef*Q_n(i-2, j)    &
                             + 10.0_rDef*Q_n(i-3, j)    &
                             -  5.0_rDef*Q_n(i-4, j)    &
                             +  1.0_rDef*Q_n(i-5, j))

            i = iMax-1
            Dis(i, j) = KCFac*( 5.0_rDef*Q_n(i+1, j)    & 
                             - 26.0_rDef*Q_n(i    , j)    &
                             + 55.0_rDef*Q_n(i-1, j)    &
                             - 60.0_rDef*Q_n(i-2, j)    &
                             + 35.0_rDef*Q_n(i-3, j)    &
                             - 10.0_rDef*Q_n(i-4, j)    &
                             +  1.0_rDef*Q_n(i-5, j))

            i = iMax-2
            Dis(i, j) = KCFac*(-10.0_rDef*Q_n(i+2, j)    & 
                             +  55.0_rDef*Q_n(i+1, j)    &
                             - 126.0_rDef*Q_n(i    , j)    &
                             + 155.0_rDef*Q_n(i-1, j)    &
                             - 110.0_rDef*Q_n(i-2, j)    &
                             +  45.0_rDef*Q_n(i-3, j)    &
                             -  10.0_rDef*Q_n(i-4, j)    &
                             +   1.0_rDef*Q_n(i-5, j))

            i = iMax-3
            Dis(i, j) = KCFac*( 10.0_rDef*Q_n(i+3, j)    & 
                             -  60.0_rDef*Q_n(i+2, j)    &
                             + 155.0_rDef*Q_n(i+1, j)    &
                             - 226.0_rDef*Q_n(i    , j)    &
                             + 205.0_rDef*Q_n(i-1, j)    &
                             - 120.0_rDef*Q_n(i-2, j)    &
                             +  45.0_rDef*Q_n(i-3, j)    &
                             -  10.0_rDef*Q_n(i-4, j)    &
                             +   1.0_rDef*Q_n(i-5, j))

            i = iMax-4
            Dis(i, j) = KCFac*( -5.0_rDef*Q_n(i+4, j)    & 
                             +  35.0_rDef*Q_n(i+3, j)    &
                             - 110.0_rDef*Q_n(i+2, j)    &
                             + 205.0_rDef*Q_n(i+1, j)    &
                             - 251.0_rDef*Q_n(i    , j)    &
                             + 210.0_rDef*Q_n(i-1, j)    &
                             - 120.0_rDef*Q_n(i-2, j)    &
                             +  45.0_rDef*Q_n(i-3, j)    &
                             -  10.0_rDef*Q_n(i-4, j)    &
                             +   1.0_rDef*Q_n(i-5, j))
       END DO
    END IF

    END SUBROUTINE Dissipation_RHS

END MODULE GetDissipation
