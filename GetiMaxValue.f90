MODULE GetiMaxValue

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: iMaxValue

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         
    SUBROUTINE iMaxValue(   xMin        ,&
                            xMax        ,&
                            x0          ,&
                            iMin        ,&
                            iMax        ,&
                            dX_Left     ,&
                            dX_Mid      ,&
                            dX_Right    ,&
                            xNew        ,&
                            deltaX_New)

    INTEGER, INTENT(IN) ::  iMin            ,&
                            xMin            ,&
                            xMax

    INTEGER, INTENT(INOUT)::    iMax    ,&
                                x0

    REAL(KIND = rDef), INTENT(INOUT) ::   dX_Left     ,&
                                        dX_Mid      ,&
                                        dX_Right

    REAL(KIND = rDef), INTENT(INOUT), DIMENSION(:), ALLOCATABLE ::    xNew    ,&
                                                                    deltaX_New

    INTEGER ::  iMin_Left   ,& 
                iMax_Left   ,&
                iMin_Mid    ,&
                iMax_Mid    ,&
                iMin_Right  ,&
                iMax_Right  ,&
                i

    x0          = 1

    dX_Left     = 0.2_rDef
    dX_Mid      = 0.006451612903226_rDef
    dX_Right    = 0.2_rDef
    
    iMin_Left   = 1
    iMax_Left   = (xMax-x0)/(dX_Right) + 1
    
    iMin_Mid    = iMax_Left
    iMax_Mid    = iMin_Mid+NINT((x0+x0)/(dX_Mid))

    iMin_Right  = iMax_Mid+1
    iMax_Right  = iMax_Mid + (-x0-xMin)/(dX_Left)

    iMax        = iMax_Right
    
    ALLOCATE(   deltaX_New(iMax)  ,&
                xNew(iMax))

    DO i = iMin, iMax
        xNew(i)         = 0.0_rDef
        deltaX_New(i)   = 0.0_rDef
    END DO

    DO i = iMin, iMax
        IF ((i >= iMin_Left).AND.(i <= iMax_Left)) THEN
            deltaX_New(i) = dX_Right 
        ELSEIF ((i > iMin_Mid).AND.(i <= iMax_Mid)) THEN 
            deltaX_New(i) = dX_Mid
        ELSEIF ((i > iMin_Right).AND.(i <= iMax_Right)) THEN 
            deltaX_New(i) = dX_Left
        END IF
    END DO

    DO i = iMin, iMax
        IF (i == iMin) THEN
            xNew(i) = REAL(xMin, rDef)
        ELSE
            xNew(i) = xNew(i-1) + deltaX_New(i)
        END IF
    END DO
    END SUBROUTINE iMaxValue

END MODULE GetiMaxValue
