MODULE GetDX

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: DX

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         
    SUBROUTINE DX(  xMin            ,&
                    xMax            ,&
                    iMin            ,&
                    iMax            ,&
                    delta_xNew      ,& 
                    xNew            ,&
                    x0              ,&
                    dX_Left         ,&
                    dX_Mid          ,&
                    dX_Right)

    INTEGER, INTENT(IN) ::  iMin        ,&
                            iMax        ,&
                            x0          ,&
                            xMin        ,&
                            xMax

    REAL(KIND = rDef), INTENT(IN) ::    dX_Left     ,&
                                        dX_Mid      ,&
                                        dX_Right

    REAL(KIND = rDef), DIMENSION(:), INTENT(INOUT) ::   delta_xNew  ,&
                                                        xNew
    
    INTEGER ::  i

    DO i = iMin, iMax
        IF (xNew(i) < REAL(-x0, rDef)) THEN
            delta_xNew(i) = 0.0767720_rDef
        ELSEIF ((xNew(i) >= REAL(-x0, rDef)).AND.(xNew(i) < 0.0_rDef)) THEN
            delta_xNew(i)   = ((TANH(-10.0_rDef*(xNew(i) + 0.50_rDef))+&
                                0.999909204262595_rDef))&
                                *0.034550753078678_rDef+0.00767720_rDef
        ELSEIF ((xNew(i) >= 0.0_rDef).AND.(xNew(i) <= REAL(x0, rDef))) THEN
            delta_xNew(i)   = ((TANH(-10.0_rDef*(-xNew(i) + 0.50_rDef))+&
                                0.999909204262595_rDef))&
                                *0.034550753078678_rDef+0.00767720_rDef
        ELSEIF ((xNew(i) > REAL(x0, rDef)).AND.(xNew(i) <= REAL(xMax, rDef))) THEN
            delta_xNew(i) = 0.0767720_rDef
        ELSE
            delta_xNew(i) = 0.0_rDef
        END IF
    END DO

    DO i = iMin, iMax
        xNew(i) = 0.0_rDef
    END DO

    DO i = iMin, iMax
        IF (i == iMin) THEN
            xNew(i) = REAL(xMin, rDef)
        ELSE
            xNew(i) = xNew(i-1) + delta_xNew(i)
        END IF
    END DO

    END SUBROUTINE DX

END MODULE GetDX
