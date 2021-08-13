MODULE GetRHSVector
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE GetdEdXVector
    USE GetDissipation
    USE GetSpatialNOdQdT

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: RHSVector

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE RHSVector(   iMin            ,&
                            iMax            ,&
                            bMin            ,&
                            bMax            ,&
                            Q_n             ,&
                            Q_np1_k         ,&
                            Q_nm1           ,&
                            dQdTau          ,&
                            RHS             ,&
                            time            ,&
                            delta_x         ,&
                            delta_Zi        ,&
                            delta_tau       ,&
                            delta_t         ,&
                            DS              ,&
                            nD              ,&
                            Dis             ,&
                            bet             ,&
                            nTau            ,&
                            nStages         ,&
                            Newtonl2Res     ,&
                            Source_Fac      ,&
                            gam             ,&
                            Mach            ,&
                            u               ,&
                            p               ,&
                            c               ,&
                            gm1             ,&
                            rho             ,&
                            rho_Static      ,&
                            P_Static        ,&
                            p0In            ,&
                            rho0In          ,&
                            Q_Exact         ,&
                            Jac_Curv        ,&
                            Inv_Jac_Curv    ,&
                            dZidX           ,&
                            dZidt           ,&
                            nFlow           ,&
                            Area)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin     ,&
                           DS       ,&
                           nD       ,&
                           nStages  ,&
                           nFlow
    
    INTEGER, INTENT(INOUT):: nTau

    REAL(KIND = rDef), INTENT(IN):: delta_x     ,&
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
                                        Newtonl2Res, &
                                        rho_Static  ,&
                                        P_Static    ,&
                                        p0In        ,&
                                        rho0In
    
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  bet             ,&
                                                    Source_Fac      ,&
                                                    Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt           ,&
                                                    Area

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT) ::    RHS         ,&
                                                            dQdTau      ,&
                                                            Dis         ,&
                                                            Q_n         ,&
                                                            Q_Exact     ,&
                                                            Q_np1_k     ,&
                                                            Q_nm1

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  dEdX            ,&
                                                        dQdT            ,&
                                                        Q_n0            ,&
                                                        E               ,&
                                                        S

    REAL(KIND = rDef) ::    l2Fac       ,&
                            DissiFac    ,&
                            Mean_p      ,&
                            Mean_rho    ,&
                            Mean_u      ,&
                            Mean_c      ,&
                            Mean_Mach   ,&
                            Contra_U    ,&
                            E_tot       ,&
                            AreaFac


    INTEGER:: i, j
    
    ALLOCATE(   E(iMax, bMax)       ,&
                dEdX(iMax, bMax)    ,&
                dQdT(iMax, bMax)    ,&
                Q_n0(iMax, bMax)    ,&
                S(iMax, bMax))

    DissiFac = 1.0_rDef/1.0_rDef
    
! Initialize E such that it is easier to calculate dEdX
    DO i = iMin, iMax
        rho     = Q_np1_k(i, 1)*Jac_Curv(i)
        u       = Q_np1_k(i, 2)/Q_np1_k(i, 1)
        p       = (gam-1.0_rDef)*(Q_np1_k(i, 3) - &
                    0.5_rDef*Q_np1_k(i, 2)*Q_np1_k(i, 2)/&
                        Q_np1_k(i, 1))*Jac_Curv(i)              
        E_tot   = Q_np1_k(i, 3)*Jac_Curv(i)
        c       = SQRT(gam*p/rho)
        Mach    = u/c
        AreaFac = Source_Fac(i)/Area(i)

        Contra_U = dZidt(i) + dZidX(i)*u

        E(i, 1) = Q_np1_k(i, 2)*dZidX(i)

        E(i, 2) = (Q_np1_k(i, 2)*Q_np1_k(i, 2)/Q_np1_k(i, 1) + (gam-1.0_rDef)*&
                    (Q_np1_k(i, 3) - Q_np1_k(i, 2)*Q_np1_k(i, 2)*0.5_rDef/&
                    Q_np1_k(i, 1)))*dZidX(i)

        E(i, 3) = (Q_np1_k(i, 2)/Q_np1_k(i, 1))*&
                    (Q_np1_k(i, 3) + (gam-1.0_rDef)*&
                    (Q_np1_k(i, 3) - Q_np1_k(i, 2)*Q_np1_k(i, 2)*0.5_rDef/&
                    Q_np1_k(i, 1)))*dZidX(i)

        S(i, 1) = AreaFac*Q_np1_k(i, 2)
        S(i, 2) = AreaFac*Q_np1_k(i, 2)*Q_np1_k(i, 2)/Q_np1_k(i, 1)
        S(i, 3) = AreaFac*(Q_np1_k(i, 2)/Q_np1_k(i, 1))*&
                    (Q_np1_k(i, 3) + (gam-1.0_rDef)*&
                    (Q_np1_k(i, 3) - Q_np1_k(i, 2)*Q_np1_k(i, 2)*0.5_rDef/&
                    Q_np1_k(i, 1)))
    END DO
    CALL Dissipation_RHS(   iMin        = iMin         ,&
                            iMax        = iMax         ,&
                            bMin        = bMin         ,&
                            bMax        = bMax         ,&
                            delta_X     = delta_X      ,&
                            delta_tau   = delta_tau    ,&
                            nD          = nD           ,&
                            Dis         = Dis          ,&
                            Q_n         = Q_np1_k) 

    !CALL SpatialNOdQdT( iMin      = iMin      ,&
    !                    iMax      = iMax      ,&
    !                    bMin      = bMin      ,&
    !                    bMax      = bMax      ,&
    !                    E         = E         ,&
    !                    dEdX      = dEdX      ,&
    !                    delta_X   = delta_zi  ,&
    !                    delta_tau = delta_tau, &
    !                    DS        = DS)
    
    CALL dEdXVector(    iMin            = iMin            ,&
                        iMax            = iMax            ,&
                        bMin            = bMin            ,&
                        bMax            = bMax            ,&
                        E               = E               ,&
                        time            = time            ,&
                        dEdX            = dEdX            ,&
                        delta_X         = delta_X         ,&
                        delta_tau       = delta_tau       ,&
                        delta_t         = delta_t         ,&
                        DS              = DS              ,&
                        Dis             = Dis             ,&
                        gam             = gam             ,&
                        Mean_Mach       = Mean_Mach       ,&
                        Mean_u          = Mean_u          ,&
                        Mean_p          = Mean_p          ,&
                        Mean_c          = Mean_c          ,&
                        gm1             = gm1             ,&
                        Mean_rho        = Mean_rho        ,&
                        rho_Static      = rho_Static      ,&
                        P_Static        = P_Static        ,&
                        p0In            = p0In            ,&
                        rho0In          = rho0In          ,&
                        Q_np1_k         = Q_np1_k         ,&
                        Q_Exact         = Q_Exact         ,&
                        Jac_Curv        = Jac_Curv        ,&
                        Inv_Jac_Curv    = Inv_Jac_Curv    ,&
                        dZidX           = dZidX           ,&
                        dZidt           = dZidt           ,&
                        delta_Zi        = delta_Zi        ,&
                        nFlow           = nFlow) 

    DO j = bMin, bMax
        DO i = iMin, iMax
            dQdT(i, j) = (  3.0_rDef*Q_np1_k(i, j)      -&
                            4.0_rDef*Q_n(i, j)          +&
                            1.0_rDef*Q_nm1(i, j))&
                            /(2.0_rDef*delta_t)
        END DO
    END DO

    DO j = bMin, bMax
        DO i = iMin, iMax
            RHS(i, j) = -delta_t*(dQdT(i, j) + dEdX(i, j) - S(i, j)) &
                            + DissiFac*Dis(i, j)
        END DO
    END DO

    !l2Fac = 0.0_rDef
    !DO j = bMin, bMax
    !    DO i = iMin, iMax
    !        l2Fac = l2Fac+RHS(i, j)*RHS(i, j)
    !    END DO
    !END DO
    !Newtonl2Res = SQRT(l2Fac)/REAL((iMax)*3, rDef)

    END SUBROUTINE RHSVector

END MODULE GetRHSVector
