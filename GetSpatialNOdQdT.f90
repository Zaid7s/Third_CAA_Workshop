MODULE GetSpatialNOdQdT

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: SpatialNOdQdT

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE SpatialNOdQdT(   iMin        ,&
                                iMax        ,&
                                bMin        ,&
                                bMax        ,&
                                E           ,&
                                dEdX        ,&
                                delta_X     ,&
                                delta_tau   ,&
                                DS           )

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS

    REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                     delta_tau
                                        
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: E

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: dEdX

    INTEGER:: i, j
    
! dEdX for Interior Domain
    DO j = bMin, bMax
        DO i = iMin, iMax 
            dEdX(i, j) = 0.0_rDef 
        END DO
    END DO

    IF      (DS == 1) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax
                IF (i == iMin) THEN
                    dEdX(i, j) = (  -3.0_rDef*E(i, j)           &
                                    +4.0_rDef*E(i+1, j)       &
                                    -         E(i+2, j))      &
                                            /(2.0_rDef*delta_x)
                ELSE IF (i == iMax) THEN
                    dEdX(i, j) = -( -3.0_rDef*E(i, j)           &
                                    +4.0_rDef*E(i-1, j)       &
                                    -         E(i-2, j))      &
                                            /(2.0_rDef*delta_x)
                ELSE
                    dEdX(i, j)   = (1.0_rDef/(delta_x*2.0_rDef))&
                                        *(E(i+1, j) - E(i-1, j))
                END IF
            END DO
        END DO
    ELSE IF (DS == 2) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax 
                IF (i == iMin) THEN
                    dEdX(i, j) = (-25.0_rDef*E(i+0, j)&
                                  +48.0_rDef*E(i+1, j)&
                                  -36.0_rDef*E(i+2, j)&
                                  +16.0_rDef*E(i+3, j)&
                                   -3.0_rDef*E(i+4, j))&
                                      /(12.0_rDef*delta_x)  
                ELSE IF (i == iMin+1) THEN
                    dEdX(i, j) =  (-3.0_rDef*E(i-1, j)&
                                  -10.0_rDef*E(i+0, j)&
                                  +18.0_rDef*E(i+1, j)&
                                   -6.0_rDef*E(i+2, j)&
                                   +1.0_rDef*E(i+3, j))&
                                    /(12.0_rDef*delta_x)  
                ELSE IF (i == iMax-1) THEN
                    dEdX(i, j) =  -(-3.0_rDef*E(i  + 1, j)&  
                                   -10.0_rDef*E(i  - 0, j)&
                                    +18.0_rDef*E(i-1, j)&
                                     -6.0_rDef*E(i-2, j)&
                                     +1.0_rDef*E(i-3, j))&
                                      /(12.0_rDef*delta_x)  
                ELSE IF (i == iMax    ) THEN
                    dEdX(i, j) = -(-25.0_rDef*E(i-0, j)&
                                   +48.0_rDef*E(i-1, j)&
                                   -36.0_rDef*E(i-2, j)&
                                   +16.0_rDef*E(i-3, j)&
                                    -3.0_rDef*E(i-4, j))&
                                      /(12.0_rDef*delta_x)  
                ELSE
                    dEdX(i, j)   = (1.0_rDef/(delta_x*12.0_rDef))&
                                            *(        E(i-2, j)&
                                            -8.0_rDef*E(i-1, j)&
                                            +8.0_rDef*E(i+1, j)&
                                            -1.0_rDef*E(i+2, j))
                END IF
            END DO
        END DO
    ELSE IF (DS == 3) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax
                IF (i == iMin) THEN
                    dEdX(i, j) = (-147.0_rDef*E(i+0, j)  &
                                 +360.0_rDef*E(i+1, j)   &
                                 -450.0_rDef*E(i+2, j)   &
                                 +400.0_rDef*E(i+3, j)   &
                                 -225.0_rDef*E(i+4, j)   &
                                  +72.0_rDef*E(i+5, j)   &
                                  -10.0_rDef*E(i+6, j))  &
                                      /(60.0_rDef*delta_x)
                ELSE IF (i == iMin+1) THEN
                    dEdX(i, j) =  (-10.0_rDef*E(i-1, j)   &
                                   -77.0_rDef*E(i+0, j)   &
                                  +150.0_rDef*E(i+1, j)   &
                                  -100.0_rDef*E(i+2, j)   &
                                   +50.0_rDef*E(i+3, j)   &
                                   -15.0_rDef*E(i+4, j)   &
                                    +2.0_rDef*E(i+5, j))  &
                                      /(60.0_rDef*delta_x)
                ELSE IF (i == iMin+2) THEN
                    dEdX(i, j) = (  2.0_rDef*E(i-2, j)   &
                                  -24.0_rDef*E(i-1, j)   &
                                  -35.0_rDef*E(i+0, j)   &
                                  +80.0_rDef*E(i+1, j)   &
                                  -30.0_rDef*E(i+2, j)   &
                                   +8.0_rDef*E(i+3, j)   &
                                   -1.0_rDef*E(i+4, j))  &
                                      /(60.0_rDef*delta_x)
                ELSE IF (i == iMax-2) THEN
                    dEdX(i, j) =  -( 2.0_rDef*E(i+2, j)    &
                                    -24.0_rDef*E(i+1, j)   &
                                    -35.0_rDef*E(i+0, j)   &
                                    +80.0_rDef*E(i-1, j)   &
                                    -30.0_rDef*E(i-2, j)   &
                                     +8.0_rDef*E(i-3, j)   &
                                     -1.0_rDef*E(i-4, j))  &
                                        /(60.0_rDef*delta_x)  
                ELSE IF (i == iMax-1) THEN
                    dEdX(i, j) =  -(-10.0_rDef*E(i+1, j)    &
                                     -77.0_rDef*E(i+0, j)   &
                                    +150.0_rDef*E(i-1, j)   &
                                    -100.0_rDef*E(i-2, j)   &
                                     +50.0_rDef*E(i-3, j)   &
                                     -15.0_rDef*E(i-4, j)   &
                                      +2.0_rDef*E(i-5, j))  &
                                         /(60.0_rDef*delta_x)
                ELSE IF (i == iMax    ) THEN
                    dEdX(i, j) = -(-147.0_rDef*E(i+0, j)   &
                                   +360.0_rDef*E(i-1, j)   &
                                   -450.0_rDef*E(i-2, j)   &
                                   +400.0_rDef*E(i-3, j)   &
                                   -225.0_rDef*E(i-4, j)   &
                                    +72.0_rDef*E(i-5, j)   &
                                    -10.0_rDef*E(i-6, j))  &
                                        /(60.0_rDef*delta_x)
                ELSE
                    dEdX(i, j)   = (1.0_rDef/(delta_x*60.0_rDef))&
                                        *(   -1.0_rDef*E(i-3, j)&
                                             +9.0_rDef*E(i-2, j)&
                                            -45.0_rDef*E(i-1, j)&
                                            +45.0_rDef*E(i+1, j)&
                                             -9.0_rDef*E(i+2, j)&
                                             +1.0_rDef*E(i+3, j))
                END IF
            END DO
        END DO
    ELSE IF (DS == 4) THEN
        DO j = bMin, bMax
            DO i = iMin, iMax
                IF (i == iMin) THEN
                    dEdX(i, j) = (-119.0_rDef*E(i+0, j)    &
                                   +296.0_rDef*E(i+1, j)   &
                                   -379.0_rDef*E(i+2, j)   &
                                   +344.0_rDef*E(i+3, j)   &
                                   -197.0_rDef*E(i+4, j)   &
                                    +64.0_rDef*E(i+5, j)   &
                                    -9.0_rDef*E(i  + 6, j))  &
                                        /(48.0_rDef*delta_x)
                ELSE IF (i == iMin+1) THEN
                    dEdX(i, j) = ( -9.0_rDef*E(i-1, j)   &
                                  -56.0_rDef*E(i-0, j)   &
                                 +107.0_rDef*E(i+1, j)   &
                                  -64.0_rDef*E(i+2, j)   &
                                  +29.0_rDef*E(i+3, j)   &
                                   -8.0_rDef*E(i+4, j)   &
                                   +1.0_rDef*E(i+5, j))  &
                                     /(48.0_rDef*delta_x)
                ELSE IF (i == iMin+2) THEN
                    dEdX(i, j) = (  1.0_rDef*E(i-2, j)    &
                                  -16.0_rDef*E(i-1, j)    &
                                  -35.0_rDef*E(i-0, j)    &
                                  +72.0_rDef*E(i+1, j)    &
                                  -29.0_rDef*E(i+2, j)    &
                                   +8.0_rDef*E(i+3, j)    &
                                   -1.0_rDef*E(i+4, j))   &
                                     /(48.0_rDef*delta_x)
               ELSE IF (i == iMax-2) THEN
                    dEdX(i, j) =  -(  1.0_rDef*E(i+2, j)    &
                                    -16.0_rDef*E(i+1, j)    &
                                    -35.0_rDef*E(i+0, j)    &
                                    +72.0_rDef*E(i-1, j)    &
                                    -29.0_rDef*E(i-2, j)    &
                                     +8.0_rDef*E(i-3, j)    &
                                     -1.0_rDef*E(i-4, j))   &
                                       /(48.0_rDef*delta_x)
                ELSE IF (i == iMax-1) THEN
                    dEdX(i, j) =  -( -9.0_rDef*E(i+1, j)   &
                                    -56.0_rDef*E(i+0, j)   &
                                   +107.0_rDef*E(i-1, j)   &
                                    -64.0_rDef*E(i-2, j)   &
                                    +29.0_rDef*E(i-3, j)   &
                                     -8.0_rDef*E(i-4, j)   &
                                     +1.0_rDef*E(i-5, j))  &
                                       /(48.0_rDef*delta_x)
                ELSE IF (i == iMax    ) THEN
                    dEdX(i, j) =  -(-119.0_rDef*E(i+0, j)   &
                                    +296.0_rDef*E(i-1, j)   &
                                    -379.0_rDef*E(i-2, j)   &
                                    +344.0_rDef*E(i-3, j)   &
                                    -197.0_rDef*E(i-4, j)   &
                                     +64.0_rDef*E(i-5, j)   &
                                     -9.0_rDef*E(i  - 6, j))  &
                                         /(48.0_rDef*delta_x)
                ELSE
                    dEdX(i, j)   = (1.0_rDef/(delta_x*48.0_rDef))&
                                        *(   -1.0_rDef*E(i-3, j)&
                                             +8.0_rDef*E(i-2, j)&
                                            -37.0_rDef*E(i-1, j)&
                                            +37.0_rDef*E(i+1, j)&
                                             -8.0_rDef*E(i+2, j)&
                                             +1.0_rDef*E(i+3, j))
                END IF
            END DO
        END DO
    END IF

    END SUBROUTINE SpatialNOdQdT

END MODULE GetSpatialNOdQdT
