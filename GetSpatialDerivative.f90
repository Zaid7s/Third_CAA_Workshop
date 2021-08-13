MODULE GetSpatialDerivative
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: SpatialDerivative

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE SpatialDerivative(   iMin            ,&
                                    iMax            ,&
                                    bMin            ,&
                                    bMax            ,&
                                    E               ,&
                                    dEdX            ,&
                                    delta_x         ,&
                                    delta_tau       ,&
                                    DS              ,&
                                    dZidX)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS

    REAL(KIND = rDef), INTENT(IN):: delta_tau  ,&
                                     delta_x

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN):: dZidX

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: E

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: dEdX

    INTEGER:: i, j
    
! dEdX for Interior Domain

    IF      (DS == 1) THEN
        DO j = bMin, bMax
            DO i = iMin+1, iMax-1 
                dEdX(i, j)   = (1.0_rDef/(delta_x*2.0_rDef))&
                                    *(E(i+1, j) - E(i-1, j))
            END DO
        END DO
    ELSE IF (DS == 2) THEN
        DO j = bMin, bMax
            DO i = iMin+2, iMax-2 
                dEdX(i, j)   = (1.0_rDef/(delta_x*12.0_rDef))&
                                    *(        E(i-2, j)&
                                    -8.0_rDef*E(i-1, j)&
                                    +8.0_rDef*E(i+1, j)&
                                    -1.0_rDef*E(i+2, j))
            END DO
        END DO
    ELSE IF (DS == 3) THEN
        DO j = bMin, bMax
            DO i = iMin+3, iMax-3
                dEdX(i, j)   = (1.0_rDef/(delta_x*60.0_rDef))&
                                    *(   -1.0_rDef*E(i-3, j)&
                                         +9.0_rDef*E(i-2, j)&
                                        -45.0_rDef*E(i-1, j)&
                                        +45.0_rDef*E(i+1, j)&
                                         -9.0_rDef*E(i+2, j)&
                                         +1.0_rDef*E(i+3, j))
            END DO
        END DO
    ELSE IF (DS == 4) THEN
        DO j = bMin, bMax
            DO i = iMin+3, iMax-3
                dEdX(i, j)   = (1.0_rDef/(delta_x*48.0_rDef))&
                                    *(   -1.0_rDef*E(i-3, j)&
                                         +8.0_rDef*E(i-2, j)&
                                        -37.0_rDef*E(i-1, j)&
                                        +37.0_rDef*E(i+1, j)&
                                         -8.0_rDef*E(i+2, j)&
                                         +1.0_rDef*E(i+3, j))
            END DO
        END DO
    END IF
    
    END SUBROUTINE SpatialDerivative

END MODULE GetSpatialDerivative
