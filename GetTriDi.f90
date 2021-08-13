MODULE GetTriDi

    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetRHSVector
    USE GetRungeKutta
    USE GetDiag
    USE BlockTriDiSolver1D

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: TriDi

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE TriDi(  iMin             ,&
                       iMax             ,&
                       bMin             ,&
                       bMax             ,&
                       Delta_Q          ,&
                       A_Jac            ,&
                       S_Jac            ,&
                       Q_n              ,&
                       Q_np1_k          ,&
                       Q_np1_kp1        ,&
                       Q_np1            ,&
                       Q_nm1            ,&
                       RHS              ,&
                       time             ,&
                       delta_x          ,&
                       delta_Zi         ,&
                       delta_tau        ,&
                       delta_t          ,&
                       DS               ,&
                       nD               ,&
                       Dis              ,&
                       bet              ,&
                       nTau             ,&
                       nStages          ,&
                       Scaling_Fac      ,&
                       E                ,&
                       E_GPT            ,&
                       Newtonl2Res      ,&
                       dQdTau           ,&
                       Source_Fac       ,&
                       gam              ,&
                       Mach             ,&
                       u                ,&
                       p                ,&
                       c                ,&
                       gm1              ,&
                       rho              ,&
                       rho_Static       ,&
                       P_Static         ,&
                       p0In             ,&
                       rho0In           ,&
                       Q_Exact          ,&
                       Jac_Curv         ,&
                       Inv_Jac_Curv     ,&
                       dZidX            ,&
                       dZidt            ,&
                       nFlow            ,&
                       Area)    

    INTEGER, INTENT(IN)::   iMin        ,&
                            iMax        ,&
                            bMax        ,&
                            bMin        ,&
                            DS          ,&
                            nD          ,&
                            nStages     ,&
                            nFlow

    INTEGER, INTENT(INOUT):: nTau

    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  bet             ,&
                                                    Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt           ,&
                                                    Area

    REAL(KIND = rDef), INTENT(IN):: time       ,&
                                    delta_x    ,&
                                    delta_Zi   ,&
                                    delta_t    ,&
                                    delta_tau  ,&
                                    Scaling_Fac
                                        
    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(IN) ::    S_Jac       ,&
                                                            A_Jac

    REAL(KIND = rDef), INTENT(INOUT):: Newtonl2Res

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS         ,&
                                                            Q_np1_k     ,&
                                                            Q_np1_kp1   ,&
                                                            Q_Exact     ,&
                                                            Dis         ,&
                                                            Q_np1       ,&
                                                            Q_nm1       ,&
                                                            Delta_Q     ,&
                                                            Q_n         ,&
                                                            E           ,&
                                                            E_GPT       ,&
                                                            dQdTau

    REAL(KIND = rDef), INTENT(INOUT)::  gam         ,&
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

    INTEGER ::   i, &
                 k, &
                 j, &
                 l, &
                 m, &
                 n

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::   RHS2

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   LowerDiag       ,&
                                                            Diag            ,&
                                                            UpperDiag

    INTEGER, DIMENSION(1):: iiMin, iiMax

    iiMin = iMin
    iiMax = iMax

    ALLOCATE(   LowerDiag(iMax, bMax, bMax), &
                Diag(iMax, bMax, bMax)      ,&
                UpperDiag(iMax, bMax, bMax), &
                RHS2(iMax, bMax))

    CALL DiagLHS( iMin        = iMin            ,& 
                  iMax        = iMax            ,& 
                  bMin        = bMin            ,& 
                  bMax        = bMax            ,& 
                  delta_Zi    = delta_Zi        ,&
                  delta_t     = delta_t         ,&
                  LowerDiag   = LowerDiag       ,& 
                  Diag        = Diag            ,& 
                  UpperDiag   = UpperDiag       ,&
                  Scaling_Fac = Scaling_Fac     ,&
                  A_Jac       = A_Jac           ,&
                  S_Jac       = S_Jac           ,&
                  Source_Fac  = Source_Fac      ,&
                  Jac_Curv    = Jac_Curv)

    CALL RungeKutta(iMin          = iMin            ,&
                    iMax          = iMax            ,&
                    bMin          = bMin            ,&
                    bMax          = bMax            ,&
                    Q_n           = Q_n             ,&
                    Q_np1_k       = Q_np1_k         ,&
                    Q_nm1         = Q_nm1           ,&
                    dQdTau        = dQdTau          ,&
                    RHS           = RHS             ,&
                    time          = time            ,&
                    delta_x       = delta_x         ,&
                    delta_Zi      = delta_Zi        ,&
                    delta_tau     = delta_tau       ,&
                    delta_t       = delta_t         ,&
                    DS            = DS              ,&
                    nD            = nD              ,&
                    Dis           = Dis             ,&
                    bet           = bet             ,&
                    nTau          = nTau            ,&
                    nStages       = nStages         ,&
                    Newtonl2Res   = Newtonl2Res     ,&
                    Source_Fac    = Source_Fac      ,&
                    gam           = gam             ,&
                    Mach          = Mach            ,&
                    u             = u               ,&
                    p             = p               ,&
                    c             = c               ,&
                    gm1           = gm1             ,&
                    rho           = rho             ,&
                    rho_Static    = rho_Static      ,&
                    P_Static      = P_Static        ,&
                    p0In          = p0In            ,&
                    rho0In        = rho0In          ,&
                    Q_Exact       = Q_Exact         ,&
                    Jac_Curv      = Jac_Curv        ,&
                    Inv_Jac_Curv  = Inv_Jac_Curv    ,&
                    dZidX         = dZidX           ,&
                    dZidt         = dZidt           ,&
                    nFlow         = nFlow           ,&
                    Area          = Area)

          
! Get the RHS. This includes the boundaries and dissipation     
    CALL RHSVector( iMin            = iMin              ,&
                    iMax            = iMax              ,&
                    bMin            = bMin              ,&
                    bMax            = bMax              ,&
                    Q_n             = Q_n               ,&
                    Q_np1_k         = Q_np1_k           ,&
                    Q_nm1           = Q_nm1             ,&
                    dQdTau          = dQdTau            ,&
                    RHS             = RHS               ,&
                    time            = time              ,&
                    delta_x         = delta_x           ,&
                    delta_Zi        = delta_Zi          ,&
                    delta_tau       = delta_tau         ,&
                    delta_t         = delta_t           ,&
                    DS              = DS                ,&
                    nD              = nD                ,&
                    Dis             = Dis               ,&
                    bet             = bet               ,&
                    nTau            = nTau              ,&
                    nStages         = nStages           ,&
                    Newtonl2Res     = Newtonl2Res       ,&
                    Source_Fac      = Source_Fac        ,&
                    gam             = gam               ,&
                    Mach            = Mach              ,&
                    u               = u                 ,&
                    p               = p                 ,&
                    c               = c                 ,&
                    gm1             = gm1               ,&
                    rho             = rho               ,&
                    rho_Static      = rho_Static        ,&
                    P_Static        = P_Static          ,&
                    p0In            = p0In              ,&
                    rho0In          = rho0In            ,&
                    Q_Exact         = Q_Exact           ,&
                    Jac_Curv        = Jac_Curv          ,&
                    Inv_Jac_Curv    = Inv_Jac_Curv      ,&
                    dZidX           = dZidX             ,&
                    dZidt           = dZidt             ,&
                    nFlow           = nFlow             ,&
                    Area            = Area)
    
    DO j = bMin, bMax
        DO i = iMin, iMax
            RHS2(i, j) = RHS(i, j)
        END DO
    END DO

!Initialize taQ and dQdTau Vectors
    
    CALL BlockTriMatrix(  numVar       = bMax,      & 
                          iSolveStart  = iiMin,     &
                          iSolveEnd    = iiMax,     &
                          aStart       = iiMin,     &
                          A            = LowerDiag, &
                          bStart       = iiMin,     &
                          B            = Diag,      &
                          cStart       = iiMin,     &
                          C            = UpperDiag, &
                          xStart       = iiMin,     &
                          X            = Delta_Q,   &
                          rStart       = iiMin,     &
                          R            = RHS)
    
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

    END SUBROUTINE TriDi
!!!!!!!!!!!!!!!!!!!!
END MODULE GetTriDi
