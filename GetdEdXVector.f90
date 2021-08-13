MODULE GetdEdXVector
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetSpatialDerivative
    USE GetInflowBC
    USE GetOutflowBC
    USE GetInflowBC_Mean_And_Pet
    USE GetOutflowBC_Mean_And_Pet
    
    IMPLICIT NONE
    PRIVATE

    PUBLIC:: dEdXVector

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         
    SUBROUTINE dEdXVector(  iMin            ,&
                            iMax            ,&
                            bMin            ,&
                            bMax            ,&
                            E               ,&
                            time            ,&
                            dEdX            ,&
                            delta_X         ,&
                            delta_tau       ,&
                            delta_t         ,&
                            DS              ,&
                            Dis             ,&
                            gam             ,&
                            Mean_Mach       ,&
                            Mean_u          ,&
                            Mean_p          ,&
                            Mean_c          ,&
                            gm1             ,&
                            Mean_rho        ,&
                            rho_Static      ,&
                            P_Static        ,&
                            p0In            ,&
                            rho0In          ,&
                            Q_np1_k         ,&
                            Q_Exact         ,&
                            Jac_Curv        ,&
                            Inv_Jac_Curv    ,&
                            dZidX           ,&
                            dZidt           ,&
                            delta_Zi        ,&
                            nFlow)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS       ,&
                           nFlow

    REAL(KIND = rDef), INTENT(INOUT):: gam         ,&
                                        Mean_Mach   ,&
                                        Mean_u      ,&
                                        Mean_p      ,&
                                        Mean_c      ,&
                                        gm1         ,&
                                        Mean_rho    ,&
                                        rho_Static  ,&
                                        P_Static    ,&
                                        p0In        ,&
                                        rho0In

    REAL(KIND = rDef), INTENT(IN):: delta_x    ,&
                                     delta_tau  ,&
                                     delta_t    ,&
                                     delta_Zi   ,&
                                     time
                                        
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::   E       ,&
                                                        Dis     ,&
                                                        Q_np1_k, &
                                                        Q_Exact

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: dEdX

    INTEGER:: i, j
    

!Inflow Boundary Using Thmpson Style BC and GPT
    DO j = bMin, bMax
        DO i = iMin, iMax 
            dEdX(i, j) = 0.0_rDef 
        END DO
    END DO

    IF (nFlow == 1) THEN
        CALL InflowBC(  iMin            = iMin             ,&
                        iMax            = iMax             ,&
                        bMin            = bMin             ,&
                        bMax            = bMax             ,&
                        E               = E                ,&
                        time            = time             ,&
                        dEdX            = dEdX             ,&
                        delta_x         = delta_Zi         ,&
                        delta_tau       = delta_tau        ,&
                        delta_t         = delta_t          ,&
                        DS              = DS               ,&
                        gam             = gam              ,&
                        Mean_Mach       = Mean_Mach        ,&
                        Mean_u          = Mean_u           ,&
                        Mean_p          = Mean_p           ,&
                        Mean_c          = Mean_c           ,&
                        gm1             = gm1              ,&
                        Mean_rho        = Mean_rho         ,&
                        rho_Static      = rho_Static       ,&
                        P_Static        = P_Static         ,&
                        p0In            = p0In             ,&
                        rho0In          = rho0In           ,&
                        Q_np1_k         = Q_np1_k          ,&
                        Q_Exact         = Q_Exact          ,&
                        Jac_Curv        = Jac_Curv         ,&
                        Inv_Jac_Curv    = Inv_Jac_Curv     ,&
                        dZidX           = dZidX)    
    ELSE IF (nFlow == 2) THEN
        CALL InflowBC_Mean_And_Pet(  iMin            = iMin             ,&
                                    iMax            = iMax             ,&
                                    bMin            = bMin             ,&
                                    bMax            = bMax             ,&
                                    E               = E                ,&
                                    time            = time             ,&
                                    dEdX            = dEdX             ,&
                                    delta_x         = delta_Zi         ,&
                                    delta_tau       = delta_tau        ,&
                                    delta_t         = delta_t          ,&
                                    DS              = DS               ,&
                                    gam             = gam              ,&
                                    Mean_Mach       = Mean_Mach        ,&
                                    Mean_u          = Mean_u           ,&
                                    Mean_p          = Mean_p           ,&
                                    Mean_c          = Mean_c           ,&
                                    gm1             = gm1              ,&
                                    Mean_rho        = Mean_rho         ,&
                                    rho_Static      = rho_Static       ,&
                                    P_Static        = P_Static         ,&
                                    p0In            = p0In             ,&
                                    rho0In          = rho0In           ,&
                                    Q_np1_k         = Q_np1_k          ,&
                                    Q_Exact         = Q_Exact          ,&
                                    Jac_Curv        = Jac_Curv         ,&
                                    Inv_Jac_Curv    = Inv_Jac_Curv     ,&
                                    dZidX           = dZidX)    
    END IF
! dEdX for Interior Domain
    CALL SpatialDerivative( iMin       =  iMin      ,&
                            iMax       =  iMax      ,&
                            bMin       =  bMin      ,&
                            bMax       =  bMax      ,&
                            E          =  E         ,&
                            dEdX       =  dEdX      ,&
                            delta_x    =  delta_Zi  ,&
                            delta_tau  =  delta_tau, &
                            DS         =  DS        ,&
                            dZidX      =  dZidX) 

!Outflow Boundary Using Thmpson Style BC and GPT
    IF (nFlow == 1) THEN
        CALL OutflowBC( iMin            = iMin             ,&
                        iMax            = iMax             ,&
                        bMin            = bMin             ,&
                        bMax            = bMax             ,&
                        E               = E                ,&
                        time            = time             ,&
                        dEdX            = dEdX             ,&
                        delta_x         = delta_Zi         ,&
                        delta_tau       = delta_tau        ,&
                        delta_t         = delta_t          ,&
                        DS              = DS               ,&
                        gam             = gam              ,&
                        Mean_Mach       = Mean_Mach        ,&
                        Mean_u          = Mean_u           ,&
                        Mean_p          = Mean_p           ,&
                        Mean_c          = Mean_c           ,&
                        gm1             = gm1              ,&
                        Mean_rho        = Mean_rho         ,&
                        rho_Static      = rho_Static       ,&
                        P_Static        = P_Static         ,&
                        p0In            = p0In             ,&
                        rho0In          = rho0In           ,&
                        Q_np1_k         = Q_np1_k          ,&
                        Q_Exact         = Q_Exact          ,&
                        Jac_Curv        = Jac_Curv         ,&
                        Inv_Jac_Curv    = Inv_Jac_Curv     ,&
                        dZidX           = dZidX)    
    ELSE IF (nFlow == 2) THEN
        CALL OutflowBC_Mean_And_Pet( iMin            = iMin             ,&
                                    iMax            = iMax             ,&
                                    bMin            = bMin             ,&
                                    bMax            = bMax             ,&
                                    E               = E                ,&
                                    time            = time             ,&
                                    dEdX            = dEdX             ,&
                                    delta_x         = delta_Zi         ,&
                                    delta_tau       = delta_tau        ,&
                                    delta_t         = delta_t          ,&
                                    DS              = DS               ,&
                                    gam             = gam              ,&
                                    Mean_Mach       = Mean_Mach        ,&
                                    Mean_u          = Mean_u           ,&
                                    Mean_p          = Mean_p           ,&
                                    Mean_c          = Mean_c           ,&
                                    gm1             = gm1              ,&
                                    Mean_rho        = Mean_rho         ,&
                                    rho_Static      = rho_Static       ,&
                                    P_Static        = P_Static         ,&
                                    p0In            = p0In             ,&
                                    rho0In          = rho0In           ,&
                                    Q_np1_k         = Q_np1_k          ,&
                                    Q_Exact         = Q_Exact          ,&
                                    Jac_Curv        = Jac_Curv         ,&
                                    Inv_Jac_Curv    = Inv_Jac_Curv     ,&
                                    dZidX           = dZidX)    
    END IF

    END SUBROUTINE dEdXVector

END MODULE GetdEdXVector
