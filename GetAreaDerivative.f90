MODULE GetAreaDerivative

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: AreaDerivative

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE AreaDerivative(  iMin      ,&
                                iMax      ,&
                                A         ,&
                                dAdX      ,&
                                delta_X   ,&
                                DS)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           DS

    REAL(KIND = rDef), INTENT(IN):: delta_x
                                        
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN):: A

    REAL(KIND = rDef), DIMENSION(:), INTENT(INOUT):: dAdX

    INTEGER:: i
    
! dAdX for Interior Domain

    IF      (DS == 1) THEN

        i = iMin
        dAdX(i) = ( -3.0_rDef*A(i)      &
                   + 4.0_rDef*A(i+1)  &
                   -          A(i+2)) &
                       /(2.0_rDef*delta_x)

        DO i = iMin+1, iMax-1 
            dAdX(i)   = (1.0_rDef/(delta_x*2.0_rDef))&
                        *(A(i+1) - A(i-1))
        END DO

        i = iMax
        dAdX(i) = -(-3.0_rDef*A(i)       &
                   + 4.0_rDef*A(i-1)   &
                   -          A(i-2))  &
                       /(2.0_rDef*delta_x)
        
    ELSE IF (DS == 2) THEN
        i       = iMin
        dAdX(i) = (-25.0_rDef*A(i-0)   &
                   +48.0_rDef*A(i+1)   &
                   -36.0_rDef*A(i+2)   &
                   +16.0_rDef*A(i+3)   &
                    -3.0_rDef*A(i+4))  &
                    /(12.0_rDef*delta_x)
        
        i       = iMin+1
        dAdX(i) = (-3.0_rDef*A(i-1)   &
                  -10.0_rDef*A(i+0)   &
                  +18.0_rDef*A(i+1)   &
                   -6.0_rDef*A(i+2)   &
                   +1.0_rDef*A(i+3))  &
                    /(12.0_rDef*delta_x)
               
        DO i = iMin+2, iMax-2 
            dAdX(i)   = (1.0_rDef/(delta_x*12.0_rDef))&
                                *(        A(i-2)&
                                -8.0_rDef*A(i-1)&
                                +8.0_rDef*A(i+1)&
                                -1.0_rDef*A(i+2))
        END DO
        
        i       = iMax-1
        dAdX(i) = -(-3.0_rDef*A(i+1)   &
                   -10.0_rDef*A(i-0)   &
                   +18.0_rDef*A(i-1)   &
                    -6.0_rDef*A(i-2)   &
                    +1.0_rDef*A(i-3))  &
                     /(12.0_rDef*delta_x)

        i       = iMax
        dAdX(i) = -(-25.0_rDef*A(i+0)   &
                    +48.0_rDef*A(i-1)   &
                    -36.0_rDef*A(i-2)   &
                    +16.0_rDef*A(i-3)   &
                     -3.0_rDef*A(i-4))  &
                     /(12.0_rDef*delta_x)

    ELSE IF (DS == 3) THEN
        
        i       = iMin
        dAdX(i) = (     -147.0_rDef*A(i+0)    &
                        +360.0_rDef*A(i+1)    &
                        -450.0_rDef*A(i+2)    &
                        +400.0_rDef*A(i+3)    &
                        -225.0_rDef*A(i+4)    &
                         +72.0_rDef*A(i+5)    &
                         -10.0_rDef*A(i+6))   &
                            /(60.0_rDef*delta_x)

        i       = iMin+1
        dAdX(i) = (     -10.0_rDef*A(i-1)     &
                        -77.0_rDef*A(i+0)     &
                       +150.0_rDef*A(i+1)     &
                       -100.0_rDef*A(i+2)     &
                        +50.0_rDef*A(i+3)     &
                        -15.0_rDef*A(i+4)     &
                         +2.0_rDef*A(i+5))    &
                            /(60.0_rDef*delta_x)
        
        i       = iMin+2
        dAdX(i) = (      2.0_rDef*A(i-2)      &
                       -24.0_rDef*A(i-1)      &
                       -35.0_rDef*A(i+0)      &
                       +80.0_rDef*A(i+1)      &    
                       -30.0_rDef*A(i+2)      &
                        +8.0_rDef*A(i+3)      &
                        -1.0_rDef*A(i+4))     &
                            /(60.0_rDef*delta_x)
        
        DO i = iMin+3, iMax-3
            dAdX(i)   = (1.0_rDef/(delta_x*60.0_rDef))&
                                *(   -1.0_rDef*A(i-3)&
                                     +9.0_rDef*A(i-2)&
                                    -45.0_rDef*A(i-1)&
                                    +45.0_rDef*A(i+1)&
                                     -9.0_rDef*A(i+2)&
                                     +1.0_rDef*A(i+3))
        END DO
        
        i       = iMax-2
        dAdX(i) = -(  2.0_rDef*A(i+2)      &
                    -24.0_rDef*A(i+1)      &
                    -35.0_rDef*A(i+0)      &
                    +80.0_rDef*A(i-1)      &    
                    -30.0_rDef*A(i-2)      &
                     +8.0_rDef*A(i-3)      &
                     -1.0_rDef*A(i-4))     &
                            /(60.0_rDef*delta_x)
        
        i       = iMax-1
        dAdX(i) = -( -10.0_rDef*A(i+1)     &
                     -77.0_rDef*A(i+0)     &
                    +150.0_rDef*A(i-1)     &
                    -100.0_rDef*A(i-2)     &
                     +50.0_rDef*A(i-3)     &
                     -15.0_rDef*A(i-4)     &
                      +2.0_rDef*A(i-5))    &
                            /(60.0_rDef*delta_x)
        
        i       = iMax
        dAdX(i) = -(-147.0_rDef*A(i-0)    &
                    +360.0_rDef*A(i-1)    &
                    -450.0_rDef*A(i-2)    &
                    +400.0_rDef*A(i-3)    &
                    -225.0_rDef*A(i-4)    &
                     +72.0_rDef*A(i-5)    &
                     -10.0_rDef*A(i-6))   &
                            /(60.0_rDef*delta_x)
    ELSE IF (DS == 4) THEN

        i = iMin
        dAdX(i) = (-119.0_rDef*A(i+0)   &
                   +296.0_rDef*A(i+1)   &
                   -379.0_rDef*A(i+2)   &
                   +344.0_rDef*A(i+3)   &
                   -197.0_rDef*A(i+4)   &
                    +64.0_rDef*A(i+5)   &
                     -9.0_rDef*A(i+6))&
                            /(48.0_rDef*delta_x)
        i = iMin+1
        dAdX(i) =  ( -9.0_rDef*A(i-1)   &
                    -56.0_rDef*A(i+0)   &
                   +107.0_rDef*A(i+1)   &
                    -64.0_rDef*A(i+2)   &
                    +29.0_rDef*A(i+3)   &
                     -8.0_rDef*A(i+4)   &
                     +1.0_rDef*A(i+5))&
                            /(48.0_rDef*delta_x)
        i = iMin+2
        dAdX(i) = (  1.0_rDef*A(i-2)   &
                   -16.0_rDef*A(i-1)   &
                   -35.0_rDef*A(i+0)   &
                   +72.0_rDef*A(i+1)   &
                   -29.0_rDef*A(i+2)   &
                    +8.0_rDef*A(i+3)   &
                    -1.0_rDef*A(i+4))&
                         /(48.0_rDef*delta_x)

        DO i = iMin+3, iMax-3
            dAdX(i)   = (1.0_rDef/(delta_x*48.0_rDef))&
                                *(   -1.0_rDef*A(i-3)&
                                     +8.0_rDef*A(i-2)&
                                    -37.0_rDef*A(i-1)&
                                    +37.0_rDef*A(i+1)&
                                     -8.0_rDef*A(i+2)&
                                     +1.0_rDef*A(i+3))
        END DO
    
        i       = iMax-2
        dAdX(i) = -( 1.0_rDef*A(i+2)   &
                   -16.0_rDef*A(i+1)   &
                   -35.0_rDef*A(i+0)   &
                   +72.0_rDef*A(i-1)   &
                   -29.0_rDef*A(i-2)   &
                    +8.0_rDef*A(i-3)   &
                    -1.0_rDef*A(i-4))&
                         /(48.0_rDef*delta_x)
        
        i       = iMax-1
        dAdX(i) =  -(-9.0_rDef*A(i+1)   &
                    -56.0_rDef*A(i+0)   &
                   +107.0_rDef*A(i-1)   &
                    -64.0_rDef*A(i-2)   &
                    +29.0_rDef*A(i-3)   &
                     -8.0_rDef*A(i-4)   &
                     +1.0_rDef*A(i-5))&
                            /(48.0_rDef*delta_x)
        
        i       = iMax
        dAdX(i) = -(-119.0_rDef*A(i-0)   &
                    +296.0_rDef*A(i-1)   &
                    -379.0_rDef*A(i-2)   &
                    +344.0_rDef*A(i-3)   &
                    -197.0_rDef*A(i-4)   &
                     +64.0_rDef*A(i-5)   &
                      -9.0_rDef*A(i-6))&
                            /(48.0_rDef*delta_x)

    END IF

    END SUBROUTINE AreaDerivative

END MODULE GetAreaDerivative
