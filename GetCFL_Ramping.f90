!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Switched Evolution Relaxation (SER) Method By Mulder and van Leer           ! 
!! SER is a CFL evolution scheme that ties CFL                                 !
!! ramping to the nonlinear residual                                           !
!!           CFL^{n+1}   = epsi*CFL^{n}                                        !
!!           epsi        = Residual^{n-1}/Residual^{n}                         !
!! There are also CFL max and min limits                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! numTimeSteps  ----> Referes to the Total Numbber of Time Steps (Time Levels)!
!! nT            ----> Referes to the current time level                       !
!! nu            ----> Is the CFL which is used both as an input and output    !
!!                       (CFL is modified after it leaves)                     !
!! Residual      ----> Referes to the Residual                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!MODULE GetCFL_Ramping
!    USE,  INTRINSIC:: ISO_FORTRAN_ENV
!    
!    IMPLICIT NONE
!    PRIVATE
!    PUBLIC:: CFL_Ramping
!
!    INTEGER, PARAMETER:: rDef = REAL64
!
!    CONTAINS
!         
!    SUBROUTINE CFL_Ramping( Residual        ,&
!                            nu              ,&
!                            numTimeSteps    ,&
!                            nT)
!
!    REAL(KIND = rDef), INTENT(INOUT)    ::  nu 
!    REAL(KIND = rDef), INTENT(IN)       :: Residual
!    INTEGER,  INTENT(IN) ::  numTimeSteps, &
!                             nT
!
!    !!! Define Local Variables !!!
!    REAL(KIND = rDef):: epsi
!    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE:: Store_Residual
!
!    ALLOCATE(Store_Residual(numTimeSteps))
!
!    Store_Residual(nT) = Residual
!
!    IF (nT >= 2) THEN
!        epsi = Store_Residual(nT-1)/Store_Residual(nT)
!    ELSE
!        epsi = 1.0_rDef
!    END IF
!
!    nu = epsi*nu
!    
!    !!!! This sets the minimum CFL !!!!!
!    IF (nu < 1.0_rDef) THEN
!        nu = 1.0_rDef
!    END IF
!    !!!! This sets the maximum CFL !!!!!
!    IF (nu > 100.0_rDef) THEN
!        nu = 100.0_rDef
!    END IF
!
!    END SUBROUTINE CFL_Ramping
!
!END MODULE GetCFL_Ramping

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Switched Evolution Relaxation (SER) Method By Mulder and van Leer           ! 
! SER is a CFL evolution scheme that ties CFL                                 !
! ramping to the nonlinear residual                                           !
!           CFL^{n+1}   = epsi*CFL^{n}                                        !
!           epsi        = Residual^{n-1}/Residual^{n}                         !
! There are also CFL max and min limits                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! numTimeSteps  ----> Referes to the Total Numbber of Time Steps (Time Levels)!
! nT            ----> Referes to the current time level                       !
! nu            ----> Is the CFL which is used both as an input and output    !
!                       (CFL is modified after it leaves)                     !
! Residual      ----> Referes to the Residual                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Test

MODULE GetCFL_Ramping
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: CFL_Ramping, &
        CFLRampObject

    INTEGER, PARAMETER:: rDef = REAL64

    INTERFACE CFL_Ramping
        MODULE PROCEDURE CFL_Ramp
    END INTERFACE CFL_Ramping

        TYPE CFLRampObject
            PRIVATE

            REAL(KIND = rDef)    ::  nu 
            REAL(KIND = rDef)    :: Residual
            INTEGER   ::  numTimeSteps, &
                nT

        END TYPE CFLRampObject

        !!! Define Local Variables !!!
        REAL(KIND = rDef):: epsi
        REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE:: Store_Residual
!!
    CONTAINS

        SUBROUTINE CFL_Ramp( object, &
        CFL)

        TYPE(CFLRampObject), INTENT(INOUT):: object

        REAL(KIND = rDef), INTENT(INOUT) ::  CFL

        object%nu = CFL

        ! SUBROUTINE CFL_Ramping( Residual        ,&
   !                         nu              ,&
   !                         numTimeSteps    ,&
   !                         nT)

   ! REAL(KIND = rDef), INTENT(INOUT)    ::  nu 
   ! REAL(KIND = rDef), INTENT(IN)       :: Residual
   ! INTEGER,  INTENT(IN) ::  numTimeSteps, &
   !                          nT

   ! !!! Define Local Variables !!!
   ! REAL(KIND = rDef):: epsi
   ! REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE:: Store_Residual

        ALLOCATE(Store_Residual(object%numTimeSteps))

        Store_Residual(object%nT) = object%Residual

        IF (object%nT >= 2) THEN
            epsi = Store_Residual(object%nT-1)/Store_Residual(object%nT)
        ELSE
            epsi = 1.0_rDef
        END IF

        object%nu = epsi*object%nu

    !!!! This sets the minimum CFL !!!!!
        IF (object%nu < 1.0_rDef) THEN
            object%nu = 1.0_rDef
        END IF
        !!!! This sets the maximum CFL !!!!!
        IF (object%nu > 100.0_rDef) THEN
            object%nu = 100.0_rDef
        END IF

    END SUBROUTINE CFL_Ramp

END MODULE GetCFL_Ramping
