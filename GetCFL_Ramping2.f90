!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Switched Evolution Relaxation (SER) Method By Mulder and van Leer
! SER is a CFL evolution scheme that ties CFL ramping to the nonlinear residual
!           CFL^{n+1}   = epsi*CFL^{n}
!           epsi        = Residual^{n-1}/Residual^{n}
! There are also CFL max and min limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! numTimeSteps  ----> Referes to the Total Numbber of Time Steps (Time Levels)
! nT            ----> Referes to the current time level
! nu            ----> Is the CFL which is used both as an input and output
!                       (CFL is modified after it leaves)
! Residual      ----> Referes to the Residual
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE GetCFL_Ramping2

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: CFL_Ramping2

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         
    SUBROUTINE CFL_Ramping2(    Residual1       ,&
                                Residual2       ,&
                                nu)

    REAL(KIND = rDef), INTENT(INOUT)    ::  nu 
    REAL(KIND = rDef), INTENT(IN)       ::  Residual1    ,&
                                            Residual2

    !!! Define Local Variables !!!
    REAL(KIND = rDef):: epsi

    epsi = Residual2/Residual1
    nu = epsi*nu
    
    !!!! This sets the minimum CFL !!!!!
    IF (nu < 1.0_rDef) THEN
        nu = 1.0_rDef
    END IF
    !!!! This sets the maximum CFL !!!!!
    IF (nu > 100.0_rDef) THEN
        nu = 100.0_rDef
    END IF


    END SUBROUTINE CFL_Ramping2

END MODULE GetCFL_Ramping2
