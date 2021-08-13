MODULE GetLUSGS
!!! Working Test Case
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetRHSVector
    USE GetRungeKutta
    USE GetLMatrix
    USE GetUMatrix
    USE GetDMatrix
    USE ForwardSweep
    USE BackwardSweep

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: LUSGS

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE LUSGS( iMin          ,&
                      iMax          ,&
                      bMin          ,&
                      bMax          ,&
                      Q_n           ,&
                      Q_np1_k       ,&
                      Q_np1         ,&
                      Q_nm1         ,&
                      Q_np1_kp1     ,&
                      RHS           ,&
                      time          ,&
                      dQdTau        ,&
                      delta_x       ,&
                      delta_Zi      ,&
                      delta_tau     ,&
                      delta_t       ,&
                      DS            ,&
                      nD            ,&
                      Dis           ,&
                      L_LowerDiag   ,&
                      L_Diag        ,&
                      U_UpperDiag   ,&
                      U_Diag        ,&
                      D_Diag        ,&
                      A_Plus        ,&
                      A_Minus       ,&
                      S_Plus        ,&
                      S_Minus       ,&
                      S_Jac         ,&
                      Delta_Q_star2, &
                      AA            ,&
                      BB            ,&
                      EE            ,&
                      Delta_Q_star  ,&
                      Delta_Q       ,&
                      bet           ,&
                      nTau          ,&
                      nStages       ,&
                      nI            ,&
                      Newtonl2Res   ,&
                      Source_Fac    ,&
                      gam           ,&
                      Mach          ,&
                      u             ,&
                      p             ,&
                      c             ,&
                      gm1           ,&
                      rho           ,&
                      rho_Static    ,&
                      P_Static      ,&
                      p0In          ,&
                      rho0In        ,&
                      Q_Exact       ,&
                      Jac_Curv      ,&
                      Inv_Jac_Curv  ,&
                      dZidX         ,&
                      dZidt         ,&
                      nFlow         ,&
                      Area)    

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS       ,&
                           nD       ,&
                           nStages  ,&
                           nI       ,&
                           nFlow
    
    INTEGER, INTENT(INOUT):: nTau

    REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                     delta_t    ,&
                                     delta_Zi   ,&
                                     delta_tau  ,&
                                     time

    REAL(KIND = rDef), INTENT(INOUT):: gam         ,&
                                        Mach        ,&
                                        u           ,&
                                        p           ,&
                                        c           ,&
                                        gm1         ,&
                                        rho         ,&
                                        rho_Static  ,&
                                        P_Static    ,&
                                        p0In        ,&
                                        rho0In

    REAL(KIND = rDef), INTENT(INOUT):: Newtonl2Res
    
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  bet             ,&
                                                    Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt           ,&
                                                    Area
    
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS             ,&
                                                            Q_np1_k         ,&
                                                            Q_n             ,&
                                                            Q_Exact         ,&
                                                            dQdTau          ,&
                                                            Q_np1           ,&
                                                            Q_nm1           ,&
                                                            Q_np1_kp1       ,&
                                                            Dis             ,&
                                                            Delta_Q_star2   ,&
                                                            Delta_Q_star    ,&
                                                            Delta_Q         ,&
                                                            AA              ,&
                                                            BB              ,&
                                                            EE

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(INOUT):: L_LowerDiag     ,&
                                                            L_Diag          ,&
                                                            U_UpperDiag     ,&
                                                            U_Diag          ,&
                                                            D_Diag          ,&
                                                            A_Minus         ,&
                                                            A_Plus          ,&
                                                            S_Jac           ,&
                                                            S_Minus         ,&
                                                            S_Plus
    
    REAL(KIND = rDef):: DissiFac, l2Fac

    REAL(KIND = rDef), DIMENSION(iMax, bMax):: dQdT        ,&
                                                dEdX        ,&
                                                E           ,&
                                                RHS2

    INTEGER            ::   i               ,&
                            k               ,&
                            j               ,&
                            l               ,&
                            m               ,&
                            n
    LOGICAL:: debug

    debug = .TRUE.
    debug = .FALSE.

! The system of equations solved here is in the form of Ax = B. Using the LUSGS 
! we break down the A Matrix to three matrices:
!                                       - L      Matrix 
!                                       - (D^-1) Matrix
!                                       - U      Matrix
! The B Vector is the RHS. The X Vector is Delta_Q. The solution for Delta_Q
! is obtained over three steps:
!               (1)  L{Delta_Q^**}    = RHS
!               (2)  D^(-1){Delta_Q*} = {Delta_Q^**}
!               (3)  U{Delta_Q}       = {Delta_Q^*}

! Step (1) is easily solved by a Forward Sweep. 
! Step (2) is solved by a direct multiplication since :
!                                       {{D^-1}^-1}{Delta_Q^**} == D{Delta_Q^**}
! Step (3) is solved by a Backward Sweep

! This subroutine creates the L Matrix 
    CALL    LMatrix( iMin         = iMin            ,&
                     iMax         = iMax            ,&
                     bMin         = bMin            ,&
                     bMax         = bMax            ,&
                     delta_x      = delta_x         ,&
                     delta_Zi     = delta_Zi        ,&
                     delta_t      = delta_t         ,&
                     A_Plus       = A_Plus          ,&
                     A_Minus      = A_Minus         ,&
                     S_Plus       = S_Plus          ,&
                     S_Minus      = S_Minus         ,&
                     S_Jac        = S_Jac           ,&
                     Source_Fac   = Source_Fac      ,&
                     L_LowerDiag  = L_LowerDiag     ,&
                     L_Diag       = L_Diag          ,&
                     Jac_Curv     = Jac_Curv        ,&
                     Inv_Jac_Curv = Inv_Jac_Curv    ,&
                     dZidX        = dZidX           ,&
                     dZidt        = dZidt)
    
! This subroutine creates the U Matrix 
    CALL    UMatrix( iMin           = iMin           ,&
                     iMax           = iMax           ,&
                     bMin           = bMin           ,&
                     bMax           = bMax           ,&
                     delta_x        = delta_x        ,&
                     delta_Zi       = delta_Zi       ,&
                     delta_t        = delta_t        ,&
                     A_Plus         = A_Plus         ,&
                     A_Minus        = A_Minus        ,&
                     S_Plus         = S_Plus         ,&
                     S_Minus        = S_Minus        ,&
                     S_Jac          = S_Jac          ,&
                     Source_Fac     = Source_Fac     ,&
                     U_UpperDiag    = U_UpperDiag    ,&
                     U_Diag         = U_Diag         ,&
                     Jac_Curv       = Jac_Curv       ,&
                     Inv_Jac_Curv   = Inv_Jac_Curv   ,&
                     dZidX          = dZidX          ,&
                     dZidt          = dZidt)

! This subroutine creates the D Matrix     
    CALL    DMatrix( iMin           = iMin           ,&
                     iMax           = iMax           ,&
                     bMin           = bMin           ,&
                     bMax           = bMax           ,&
                     delta_x        = delta_x        ,&
                     delta_Zi       = delta_Zi       ,&
                     delta_t        = delta_t        ,&
                     A_Plus         = A_Plus         ,&
                     A_Minus        = A_Minus        ,&
                     S_Plus         = S_Plus         ,&
                     S_Minus        = S_Minus        ,&
                     S_Jac          = S_Jac          ,&
                     Source_Fac     = Source_Fac     ,&
                     D_Diag         = D_Diag         ,&
                     Jac_Curv       = Jac_Curv       ,&
                     Inv_Jac_Curv   = Inv_Jac_Curv   ,&
                     dZidX          = dZidX          ,&
                     dZidt          = dZidt) 

    IF (debug) WRITE(0, *) "Survive Before RK"

    CALL RungeKutta(iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    Q_n             = Q_n           ,&
                    Q_np1_k         = Q_np1_k       ,&
                    Q_nm1           = Q_nm1         ,&
                    dQdTau          = dQdTau        ,&
                    RHS             = RHS           ,&
                    time            = time          ,&
                    delta_x         = delta_x       ,&
                    delta_Zi        = delta_Zi      ,&
                    delta_tau       = delta_tau     ,&
                    delta_t         = delta_t       ,&
                    DS              = DS            ,&
                    nD              = nD            ,&
                    Dis             = Dis           ,&
                    bet             = bet           ,&
                    nTau            = nTau          ,&
                    nStages         = nStages       ,&
                    Newtonl2Res     = Newtonl2Res   ,&
                    Source_Fac      = Source_Fac    ,&
                    gam             = gam           ,&
                    Mach            = Mach          ,&
                    u               = u             ,&
                    p               = p             ,&
                    c               = c             ,&
                    gm1             = gm1           ,&
                    rho             = rho           ,&
                    rho_Static      = rho_Static    ,&
                    P_Static        = P_Static      ,&
                    p0In            = p0In          ,&
                    rho0In          = rho0In        ,&
                    Q_Exact         = Q_Exact       ,&
                    Jac_Curv        = Jac_Curv      ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv  ,&
                    dZidX           = dZidX         ,&
                    dZidt           = dZidt         ,&
                    nFlow           = nFlow         ,&
                    Area            = Area)

        IF (debug) WRITE(0, *) "Survive After RK"
    
! Get the RHS. This includes the boundaries and dissipation     


    CALL RHSVector( iMin            = iMin          ,&
                    iMax            = iMax          ,&
                    bMin            = bMin          ,&
                    bMax            = bMax          ,&
                    Q_n             = Q_n           ,&
                    Q_np1_k         = Q_np1_k       ,&
                    Q_nm1           = Q_nm1         ,&
                    dQdTau          = dQdTau        ,&
                    RHS             = RHS           ,&
                    time            = time          ,&
                    delta_x         = delta_x       ,&
                    delta_Zi        = delta_Zi      ,&
                    delta_tau       = delta_tau     ,&
                    delta_t         = delta_t       ,&
                    DS              = DS            ,&
                    nD              = nD            ,&
                    Dis             = Dis           ,&
                    bet             = bet           ,&
                    nTau            = nTau          ,&
                    nStages         = nStages       ,&
                    Newtonl2Res     = Newtonl2Res   ,&
                    Source_Fac      = Source_Fac    ,&
                    gam             = gam           ,&
                    Mach            = Mach          ,&
                    u               = u             ,&
                    p               = p             ,&
                    c               = c             ,&
                    gm1             = gm1           ,&
                    rho             = rho           ,&
                    rho_Static      = rho_Static    ,&
                    P_Static        = P_Static      ,&
                    p0In            = p0In          ,&
                    rho0In          = rho0In        ,&
                    Q_Exact         = Q_Exact       ,&
                    Jac_Curv        = Jac_Curv      ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv  ,&
                    dZidX           = dZidX         ,&
                    dZidt           = dZidt         ,&
                    nFlow           = nFlow         ,&
                    Area            = Area)
    IF (debug) WRITE(0, *) "Survive After RHS"

    DO j = bMin, bMax
        DO i = iMin, iMax
            RHS2(i, j) = RHS(i, j)
        END DO
    END DO

! Start solving with a Forward Sweep, Step (1) 
    CALL F_Sweep (L_LowerDiag   = L_LowerDiag   ,&
                  L_Diag        = L_Diag        ,&
                  iMin          = iMin          ,&
                  iMax          = iMax          ,&
                  bMin          = bMin          ,&
                  bMax          = bMax          ,&
                  RHS           = RHS           ,&
                  Delta_Q_star2 = Delta_Q_star2)
        
! Solve the D Matrix, No Sweeps Needed Step(2)
    DO i = iMin, iMax
        DO m = iMin, bMax
            DO n = iMin, bMin
                DO l = iMin, bMax
                    AA(m, l)     = D_Diag(m, l, i)
                    BB(l, n)     = Delta_Q_star2(i, l)
                END DO 
            END DO
        END DO
        EE = MATMUL(AA, BB)
        Delta_Q_star(i, 1) = EE(1, 1)
        Delta_Q_star(i, 2) = EE(2, 1)
        Delta_Q_star(i, 3) = EE(3, 1)
    END DO
     
! Start solving with a Backward Sweep, Step (3) 
    CALL B_Sweep (U_UpperDiag   = U_UpperDiag   ,&
                  U_Diag        = U_Diag        ,&
                  iMin          = iMin          ,&
                  iMax          = iMax          ,&
                  bMin          = bMin          ,&
                  bMax          = bMax          ,&
                  Delta_Q_star  = Delta_Q_star  ,&
                  Delta_Q       = Delta_Q) 

!! Update the Solution: Q_i^{n+1} = \Delta{Q_{i}} + Q_{i} 
!!                      Q_i^{n}   = Q_i^{n+1}

    
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_np1_kp1(i, j) = Delta_Q(i, j) + Q_np1_k(i, j)
        END DO
    END DO

    DO j = bMin, bMax
        DO i = iMin, iMax
            RHS(i, j) = RHS2(i, j)
        END DO
    END DO
    
    END SUBROUTINE LUSGS

END MODULE GetLUSGS
