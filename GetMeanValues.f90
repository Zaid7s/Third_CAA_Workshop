MODULE GetMeanValues
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetLUSGS
    USE GetJacobians
    USE GetTriDi
    USE GetAreaDerivative
    USE GetDX
    USE GetCFL_Ramping
    USE GetCFL_Ramping2
    USE GetExactSol
    USE GetiMaxValue

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: MeanValues

    TYPE(CFLRampObject) :: thisCFLObject

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE MeanValues(  Q_Mean      ,&
                            Q_Exact     ,&
                            nFlow       ,&
                            p_Exit      ,&
                            nL          ,&
                            nD          ,&
                            DS)

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE, INTENT(INOUT) ::        &
                                                        Q_Mean              ,&
                                                        Q_Exact

    INTEGER, INTENT(IN):: nFlow

    INTEGER, INTENT(OUT)::  nL  ,&
                            nD  ,&
                            DS

    REAL(KIND = rDef), INTENT(INOUT):: p_Exit

    INTEGER     ::      xMin            ,&
                        xMax            ,&
                        iMax3           ,&
                        iMax4           ,&
                        iMin            ,&
                        iMax_Original   ,&
                        iMax            ,&
                        nPts            ,&
                        bMin            ,&
                        bMax            ,&
                        i               ,&
                        ii              ,&
                        nT              ,&
                        nI              ,&
                        numTimeSteps    ,&
                        MaxIteration    ,&
                        k               ,&
                        j               ,&
                        l               ,&
                        m               ,&
                        n               ,&
                        nTau            ,&
                        nStages         ,&
                        nLHS            ,&
                        iQuarter        ,&
                        xC

    REAL(KIND = rDef) ::    delta_x     ,&
                            delta_x1    ,&
                            delta_x2    ,&
                            delta_x3    ,&
                            pi          ,&
                            delta_t     ,&
                            delta_Zi    ,&
                            delta_tau   ,&
                            tol         ,&
                            tol2        ,&
                            gam         ,&
                            gm1         ,&
                            omega1      ,&
                            omega2      ,&
                            time        ,&
                            L2Fac       ,&
                            NewtonL2Fac, &
                            epsi_A      ,&
                            epsi_S      ,&
                            Scaling_Fac, &
                            kDeltaX     ,&
                            sigma_Inv   ,&
                            nu_physical, &
                            nu          ,&
                            xx          ,&
                            rho         ,&
                            u           ,&
                            p           ,&
                            Mach        ,&
                            c           ,&
                            eig1_A      ,&
                            eig2_A      ,&
                            eig3_A      ,&
                            eig1_S      ,&
                            eig2_S      ,&
                            eig3_S      ,&
                            T0In        ,&
                            p0In        ,&
                            rho0In      ,&
                            P_Static    ,&
                            rho_Static

                    
    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE ::     x               ,&
                                                        x1              ,&
                                                        delta_xNew      ,&
                                                        Area1           ,&
                                                        Jac_Curv        ,&
                                                        Inv_Jac_Curv    ,&
                                                        l2Res           ,&
                                                        Residual        ,&
                                                        MaxWaveSpeed    ,&
                                                        NewtonResidual  ,&
                                                        Newtonl2Res     ,&
                                                        bet             ,&
                                                        Area            ,&
                                                        dAdZi            ,&
                                                        Source_Fac      ,&
                                                        dXdt            ,&
                                                        dZidX           ,&
                                                        dZidt

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_peturbation   ,&
                                                        Q_n             ,&
                                                        Q_np1           ,&
                                                        Q_nm1           ,&
                                                        dQdTau          ,&
                                                        RHS             ,&
                                                        E_GPT           ,&
                                                        E               ,&
                                                        Q_np1_EGPT      ,&
                                                        Q_EGPT          ,&
                                                        Delta_Q_EGPT    ,&
                                                        Q_np1_k         ,&
                                                        Q_np1_kp1       ,&
                                                        Delta_Q_star2   ,&
                                                        Delta_Q_star    ,&
                                                        Delta_Q         ,&
                                                        AA              ,&
                                                        BB              ,&
                                                        EE              ,&
                                                        Dis             ,&
                                                        Exact

    REAL(KIND = rDef), DIMENSION(:, :, :), ALLOCATABLE ::   A_Jac           ,&
                                                            S_Jac           ,&
                                                            L_LowerDiag     ,&
                                                            L_Diag          ,&
                                                            U_UpperDiag     ,&
                                                            U_Diag          ,&
                                                            D_Diag          ,&
                                                            A_Plus          ,&
                                                            A_Minus         ,&
                                                            S_Plus          ,&
                                                            S_Minus         ,&
                                                            Q_Store

                                                            
    CHARACTER(LEN = 8)  :: DifferStencil
    CHARACTER(LEN = 4)  :: Dissipation
    CHARACTER(LEN = 5)  :: LeftHandSide
    CHARACTER(LEN = 25):: filename1
    CHARACTER(LEN = 32):: filename2
    CHARACTER(LEN = 38):: filename3
    CHARACTER(LEN = 17):: filename4
    CHARACTER(LEN = 12):: filename5
    CHARACTER(LEN = 18):: filename6
    CHARACTER(LEN = 18):: filename7
    
    pi              = 4.0_rDef*ATAN(1.0_rDef)
    xMin            = -10                      
    xMax            =  10                    
    delta_Zi        = 1.0_rDef     
    iMin            = 1
    bMin            = 1
    bMax            = 3
    tol             = 1E-17_rDef
    tol2            = 1E-15_rDef
    gam             = 1.40_rDef
    gm1             = 1.0_rDef/(gam-1.0_rDef)
    numTimeSteps    = 200001
    numTimeSteps    = 101
    MaxIteration    = 151
    omega1          = pi/10.0_rDef
    omega2          = pi/50.0_rDef
    time            = 0.0_rDef
    nStages         = 4
    delta_t         = 0.0_rDef

    p_Exit          = 1.0_rDef/gam  ! 0.650_rDef
    p0In            = 0.7975371184_rDef 
    rho0In          = 1.081930199_rDef 
    t0In            = 1.032_rDef 
    
    CALL iMaxValue( xMin        = xMin      ,&
                    xMax        = xMax      ,&
                    x0          = xC        ,&
                    iMin        = iMin      ,&
                    iMax        = iMax      ,&
                    dX_Left     = delta_x1  ,&
                    dX_Mid      = delta_x2  ,&
                    dX_Right    = delta_x3  ,&
                    xNew        = x         ,&
                    deltaX_New  = delta_xNew)

    nPts            = iMax

    ALLOCATE(   Area(iMax)                          ,&
                Area1(iMax)                         ,&
                x1(iMax)                            ,&
                Source_Fac(iMax)                    ,&
                Q_peturbation(iMax, bMax)           ,&
                Q_Mean(iMax, bMax)                  ,&
                Q_n(iMax, bMax)                     ,&
                Q_Store(numTimeSteps, iMax, bMax)   ,&
                Q_np1(iMax, bMax)                   ,&
                Q_nm1(iMax, bMax)                   ,&
                dQdTau(iMax, bMax)                  ,&
                dAdZi(iMax)                         ,&
                RHS(iMax, bMax)                     ,&
                l2Res(numTimeSteps)                 ,&
                Residual(numTimeSteps)              ,&
                MaxWaveSpeed(numTimeSteps)          ,&
                Newtonl2Res(MaxIteration)           ,&
                NewtonResidual(MaxIteration)        ,&
                E_GPT(iMax, bMax)                   ,&
                E(iMax, bMax)                       ,&
                Dis(iMax, bMax)                     ,&
                Q_np1_k(iMax, bMax)                 ,&
                Q_Exact(iMax, bMax)                 ,&
                Q_np1_kp1(iMax, bMax)               ,&
                L_LowerDiag(bMax, bMax, iMax)       ,&
                L_Diag(bMax, bMax, iMax)            ,&
                U_UpperDiag(bMax, bMax, iMax)       ,&
                U_Diag(bMax, bMax, iMax)            ,&
                D_Diag(bMax, bMax, iMax)            ,&
                Exact(iMax, bMax)                   ,&
                Delta_Q_star2(iMax, bMax)           ,&
                Delta_Q_star(iMax, bMax)            ,&
                Delta_Q(iMax, bMax)                 ,&
                AA(bMax, bMax)                      ,&
                BB(bMax, bMin)                      ,&
                EE(bMax, bMin)                      ,&
                bet(nStages)                        ,&
                A_Jac(bMax, bMax, iMax)             ,&
                A_Plus(bMax, bMax, iMax)            ,&
                A_Minus(bMax, bMax, iMax)           ,&
                S_Jac(bMax, bMax, iMax)             ,&
                S_Plus(bMax, bMax, iMax)            ,&
                S_Minus(bMax, bMax, iMax)           ,&
                Jac_Curv(iMax)                      ,&
                Inv_Jac_Curv(iMax)                  ,&
                dXdt(iMax)                          ,&
                dZidX(iMax)                         ,&
                dZidt(iMax))
    
    DO i = iMin, iMax
        DO j = bMin, bMax 
            DO k = bMin, bMax 
                Area(i)                         = 0.0_rDef
                Area1(i)                        = 0.0_rDef
                x1(i)                           = 0.0_rDef
                Source_Fac(i)                   = 0.0_rDef
                Q_n(i, j)                       = 0.0_rDef
                Q_np1(i, j)                     = 0.0_rDef
                Q_nm1(i, j)                     = 0.0_rDef
                dQdTau(i, j)                    = 0.0_rDef
                dAdZi(i)                        = 0.0_rDef
                RHS(i, j)                       = 0.0_rDef
                E_GPT(i, j)                     = 0.0_rDef
                E(i, j)                         = 0.0_rDef
                Dis(i, j)                       = 0.0_rDef
                Q_np1_k(i, j)                   = 0.0_rDef
                Q_np1_kp1(i, j)                 = 0.0_rDef   
                L_LowerDiag(k, j, i)            = 0.0_rDef
                L_Diag(k, j, i)                 = 0.0_rDef
                U_UpperDiag(k, j, i)            = 0.0_rDef
                U_Diag(k, j, i)                 = 0.0_rDef
                D_Diag(k, j, i)                 = 0.0_rDef
                Exact(i, j)                     = 0.0_rDef   
                Delta_Q_star2(i, j)             = 0.0_rDef   
                Delta_Q_star(i, j)              = 0.0_rDef   
                Delta_Q(i, j)                   = 0.0_rDef
                AA(k, j)                        = 0.0_rDef
                BB(j, 1)                        = 0.0_rDef
                EE(j, 1)                        = 0.0_rDef
                bet(nStages)                    = 0.0_rDef   
                A_Jac(k, j, i)                  = 0.0_rDef
                A_Plus(k, j, i)                 = 0.0_rDef
                A_Minus(k, j, i)                = 0.0_rDef
                S_Jac(k, j, i)                  = 0.0_rDef
                S_Plus(k, j, i)                 = 0.0_rDef
                S_Minus(k, j, i)                = 0.0_rDef
                Jac_Curv(i)                     = 0.0_rDef
                Inv_Jac_Curv(i)                 = 0.0_rDef
                dXdt(i)                         = 0.0_rDef
                dZidX(i)                        = 0.0_rDef
                dZidt(i)                        = 0.0_rDef
            END DO
        END DO
    END DO

    WRITE(filename5, '(a12)') 'Jacobian.dat'
    
    OPEN (64, FILE = filename5, FORM = 'FORMATTED')

    CALL DX(  xMin          = xMin          ,&
              xMax          = xMax          ,&
              iMin          = iMin          ,&
              iMax          = iMax          ,&
              delta_xNew    = delta_xNew    ,&
              xNew          = x             ,&
              x0            = xC            ,&
              dX_Left       = delta_x1      ,&
              dX_Mid        = delta_x2      ,&
              dX_Right      = delta_x3)

    DO i = iMin, iMax
        WRITE(64, *) i, delta_xNew(i)
    END DO

    CLOSE(64)
    
!! RK Coefficients
    bet(1)              = (1.0_rDef)/(4.0_rDef) 
    bet(2)              = (1.0_rDef)/(3.0_rDef) 
    bet(3)              = (1.0_rDef)/(2.0_rDef) 
    bet(4)              =  1.0_rDef 
    
    DO nL = 1, 1
!! Different LHS. 1--> Uses Lower upper Symmetric Gauss Seidel Method by 
!!                      Yoon and Jameson
!!                2--> Uses ADI on the LHS
        IF (nL == 1) THEN
            LeftHandSide = "LUSGS"
        ELSE IF (nL == 2) THEN
            LeftHandSide = "TriDi"
        ELSE
            CYCLE
        END IF
    !! Different RHS Dissipation.   1--> Second Order Dissipation
    !!                              2--> Fourth Order
    !!                              3--> Sixth Order
    !!                              4--> 8th Order
    !!                              5--> Tenth Order
        DO nD = 5, 5
            IF (nD == 1) THEN
                Dissipation = 'D_02'
            ELSE IF (nD == 2) THEN
                Dissipation = 'D_04'
            ELSE IF (nD == 3) THEN
                Dissipation = 'D_06'
            ELSE IF (nD == 4) THEN
                Dissipation = 'D_08'
            ELSE IF (nD == 5) THEN
                Dissipation = 'D_10'
            ELSE
                EXIT
            END IF
        !! Different RHS Differncing Stencil.   1--> Second Order
        !!                                      2--> Fourth Order
        !!                                      3--> Sixth Order
        !!                                      4--> RDRP
            DO DS = 4, 4
                IF (DS == 1) THEN
                    Scaling_Fac     = 1.8_rDef 
                    DifferStencil   = 'SecOrder'
                    kDeltaX         = 1.0_rDef
                ELSE IF (DS == 2) THEN
                    Scaling_Fac     = 3.791411_rDef
                    DifferStencil   = 'FouOrder'
                    kDeltaX         = 1.400_rDef
                ELSE IF (DS == 3) THEN
                    Scaling_Fac     = 10.567013_rDef
                    DifferStencil   = 'SixOrder'
                    kDeltaX         = 1.586_rDef
                ELSE IF (DS == 4) THEN
                    Scaling_Fac     = 10.6597009_rDef
                    DifferStencil   = 'RDRPSten'
                    kDeltaX         = 1.664_rDef
                ELSE
                    EXIT
                END IF

                IF ((DS == 1).AND.(nD < 2)) THEN
                    CYCLE
                END IF

                IF ((DS == 2).AND.(nD < 4)) THEN
                    CYCLE
                END IF

                IF ((DS == 3).AND.(nD < 5)) THEN
                    CYCLE
                END IF
                
                IF ((DS == 4).AND.(nD < 5)) THEN
                    CYCLE
                END IF
                
                WRITE(filename2, '(a5, a1, a8, a1, a4, a13)') &
                        LeftHandSide, '_', DifferStencil, '_', &
                        Dissipation, '_Converge.dat'
                
                OPEN (56, FILE = filename2, FORM = 'FORMATTED')

            !! Define Area and Domain

                DO i = iMin, iMax
                    IF (x(i) >= 0.0_rDef) THEN
                        Area(i) = 0.536572_rDef-0.198086_rDef*&
                                  EXP(-LOG(2.0_rDef)*((x(i)/0.6_rDef)**2.0_rDef))
                    ELSEIF (x(i) < 0.0_rDef) THEN
                        Area(i) = 1.0_rDef-0.661514_rDef*&
                                  EXP(-LOG(2.0_rDef)*((x(i)/0.6_rDef)**2.0_rDef))
                    END IF
                    !IF  (x(i) < -5.0_rDef) THEN
                    !    Area(i) = 0.53672_rDef
                    !ELSE IF ((x(i) >= -5.0_rDef).AND.(x(i) <= 5.0_rDef)) THEN
                    !    xx       = x(i)/5.0_rDef
                    !    Area1(i) =  ((((((-0.34746_rDef*xx&
                    !                     -0.436443_rDef)*xx)&
                    !                     +0.579100_rDef)*xx)&
                    !                     +0.872886_rDef)*xx*xx)&
                    !                     +0.331917_rDef
                    !    Area(i)  = Area1(i)
                    ! ELSE
                    !    Area(i) = 1.0_rDef
                    !END IF
                END DO

                DO i = iMin, iMax
                    dZidX(i)        = 1.0_rDef/delta_xNew(i)
                    Jac_Curv(i)     = dZidX(i) 
                    Inv_Jac_Curv(i) = 1.0_rDef/(Jac_Curv(i))
                    dXdt(i)         = 0.0_rDef 
                    dZidt(i)        = Jac_Curv(i)*dXdt(i)
                END DO

            !! Get dAdZi
                CALL AreaDerivative(    iMin    =   iMin        ,&
                                        iMax    =   iMax        ,&
                                        A       =   Area        ,&
                                        dAdX    =   dAdZi       ,&
                                        delta_x =   delta_Zi    ,&
                                        DS      =   DS)

            !! Write Out The area to get Exact Solution
                WRITE(filename1, '(a5, a8, a5, i3.3, a4)') &
                        "Area_", DifferStencil, '_NPTS', nPts, '.dat'
                OPEN (59, FILE = filename1, FORM = 'FORMATTED')

                !WRITE(59, *) nPts
                DO i = iMin, iMax
                    Source_Fac(i)  = -Jac_Curv(i)*(dAdZi(i))
                    WRITE(59, *) i, x(i), Area(i), Jac_Curv(i), delta_xNew(i)
                END DO
                CLOSE(59)

                CALL ExactSol(  iMin    = iMin          ,&
                                iMax    = iMax          ,&
                                p0In    = p0In          ,&
                                t0In    = T0In          ,&
                                p_Exit  = p_Exit        ,&
                                tol     = tol           ,&
                                ii      = i             ,&
                                x       = x             ,&
                                A       = Area          ,&
                                p       = Exact(:, 1)   ,&
                                rho     = Exact(:, 2)   ,&
                                u       = Exact(:, 3))

            !!! NOTE THE EXACT SOLUTION FILE WRITES OUT IN THE FOLLOWING
            !!! FORMAT : X(i), PRESSURE(i), DENSITY(i), VELOCITY(i). 
            !!! Why?  WHO THE FUCK KNOWS!!
            !!! Hence, Q_Exact has a weird format. 
                DO i = iMin, iMax
                    Q_Exact(i, 1) = Exact(i, 2)
                    Q_Exact(i, 2) = Exact(i, 2)*Exact(i, 3)
                    Q_Exact(i, 3) = Exact(i, 1)/(gam-1.0_rDef) +  &
                                    Exact(i, 2)*Exact(i, 3)*&
                                    Exact(i, 3)*0.5_rDef
                END DO

                WRITE(filename6, '(a6, a8, a4)') &
                        "Exact_", DifferStencil, '.dat'
                OPEN (67, FILE = filename6, FORM = 'FORMATTED')

                !WRITE(59, *) nPts
                DO i = iMin, iMax
                    WRITE(67, *) x(i), Q_Exact(i, 1)                        ,&
                                       !!!!!!!!!!!!Exact Velocity!!!!!!!!
                                       Q_Exact(i, 2)/Q_Exact(i, 1)          ,&
                                       !!!!!!!!!!!!Exact Pressure!!!!!!!!
                                       ((gam-1.0_rDef)*(Q_Exact(i, 3) &
                                       - 0.5_rDef*Q_Exact(i, 2)*&
                                       Q_Exact(i, 2)/Q_Exact(i, 1)))        ,&
                                       !!!!!!!!!!!!Exact Mach Number!!!!!
                                       (Q_Exact(i, 2)/Q_Exact(i, 1))/&
                                       (SQRT(gam*(((gam-1.0_rDef)*&
                                       (Q_Exact(i, 3) - 0.5_rDef*&
                                       Q_Exact(i, 2)*Q_Exact(i, 2)&
                                       /Q_Exact(i, 1))))/(Q_Exact(i, 1))))

                END DO
                CLOSE(67)

            !! Initial Condition
                DO i = iMin, iMax
                    rho         = 1.0_rDef  ! Exact(i, 2)
                    u           = 0.4_rDef  ! Exact(i, 3)
                    p           = p_Exit  ! Exact(i, 1)
                    Q_n(i, 1)   = rho*Inv_Jac_Curv(i)
                    Q_n(i, 2)   = rho*u*Inv_Jac_Curv(i)
                    Q_n(i, 3)   = (p/(gam-1.0_rDef))*Inv_Jac_Curv(i) +    &
                                    0.5_rDef*rho*u*u*Inv_Jac_Curv(i)
                END DO 

                DO j = bMin, bMax
                    DO i = iMin, iMax
                        Q_Exact(i, j)   = Q_Exact(i, j)*Inv_Jac_Curv(i)
                    END DO
                END DO

            !! Needed For Outflow Boundary Calculation
                P_Static    = p_Exit*Inv_Jac_Curv(1)
                rho_Static  = Q_n(1, 1)
                nu_physical = 2.5_rDef


                DO nT = 1, numTimeSteps

                    WRITE(filename3, '(a5, a1, a8, a1, a4, a1, i8.8, a2, a8)')      &
                            LeftHandSide, '_', DifferStencil, '_',           &
                            Dissipation, '_', nT, 'nT', 'Mean.dat' 

                    DO j = bMin, bMax
                        DO i = iMin, iMax
                            Q_Store(nT, i, j) = Q_n(i, j)
                        END DO
                    END DO

                    DO j = bMin, bMax
                        DO i = iMin, iMax
                            Q_np1_k(i, j)   = Q_Store(nT, i, j)
                        END DO
                    END DO

                    IF (nT == 1) THEN
                        DO j = bMin, bMax
                            DO i = iMin, iMax
                                Q_nm1(i, j)   = Q_Store(nT, i, j)
                            END DO
                        END DO
                    ELSE
                        DO j = bMin, bMax
                            DO i = iMin, iMax
                                Q_nm1(i, j)   = Q_Store(nT-1, i, j)
                            END DO
                        END DO
                    END IF

                    DO nI = 1, MaxIteration
                        CALL Jacobians(  iMin           = iMin          ,&
                                         iMax           = iMax          ,&
                                         bMin           = bMin          ,&
                                         bMax           = bMax          ,&
                                         A_Jac          = A_Jac         ,&
                                         u              = u             ,&
                                         p              = p             ,&
                                         c              = c             ,&
                                         rho            = rho           ,&
                                         Mach           = Mach          ,&
                                         Q_n            = Q_np1_k       ,&
                                         gm1            = gm1           ,&
                                         gam            = gam           ,&
                                         eig1_A         = eig1_A        ,&
                                         eig2_A         = eig2_A        ,&
                                         eig3_A         = eig3_A        ,&
                                         epsi_A         = epsi_A        ,&
                                         Source_Fac     = Source_Fac    ,&
                                         A_Plus         = A_Plus        ,&
                                         A_Minus        = A_Minus       ,&
                                         eig1_S         = eig1_S        ,&
                                         eig2_S         = eig2_S        ,&
                                         eig3_S         = eig3_S        ,&
                                         epsi_S         = epsi_S        ,&
                                         S_Plus         = S_Plus        ,&
                                         S_Minus        = S_Minus       ,&
                                         S_Jac          = S_Jac         ,&
                                         Jac_Curv       = Jac_Curv      ,&
                                         Inv_Jac_Curv   = Inv_Jac_Curv  ,&
                                         dZidX          = dZidX         ,&
                                         dZidt          = dZidt         ,&
                                         Area           = Area)          

                        sigma_Inv       = epsi_A*kDeltaX
                        nu              = 2.3_rDef
                        
                        IF (MOD(nt, 2) == 1) THEN
                            delta_t         = nu_physical/SQRT(sigma_Inv&
                                                              *sigma_Inv)
                            delta_tau       = (nu)/SQRT(sigma_Inv*sigma_Inv &
                                                + (1.0_rDef/delta_t)&
                                                *(1.0_rDef/delta_t))
                        END IF

                        time = time+delta_t

                        IF (nL == 1) THEN
                            CALL LUSGS( iMin            = iMin            ,&
                                        iMax            = iMax            ,&
                                        bMin            = bMin            ,&
                                        bMax            = bMax            ,&
                                        Q_n             = Q_n             ,&
                                        Q_np1_k         = Q_np1_k         ,&
                                        Q_np1           = Q_np1           ,&
                                        Q_nm1           = Q_nm1           ,&
                                        Q_np1_kp1       = Q_np1_kp1       ,&
                                        RHS             = RHS             ,&
                                        time            = time            ,&
                                        dQdTau          = dQdTau          ,&
                                        delta_x         = delta_x         ,&
                                        delta_Zi        = delta_Zi        ,&
                                        delta_tau       = delta_tau       ,&
                                        delta_t         = delta_t         ,&
                                        DS              = DS              ,&
                                        nD              = nD              ,&
                                        Dis             = Dis             ,&
                                        L_LowerDiag     = L_LowerDiag     ,&
                                        L_Diag          = L_Diag          ,&
                                        U_UpperDiag     = U_UpperDiag     ,&
                                        U_Diag          = U_Diag          ,&
                                        D_Diag          = D_Diag          ,&
                                        A_Plus          = A_Plus          ,&
                                        A_Minus         = A_Minus         ,&
                                        S_Plus          = S_Plus          ,&
                                        S_Minus         = S_Minus         ,&
                                        S_Jac           = S_Jac           ,&
                                        Delta_Q_star2   = Delta_Q_star2   ,&
                                        AA              = AA              ,&
                                        BB              = BB              ,&
                                        EE              = EE              ,&
                                        Delta_Q_star    = Delta_Q_star    ,&
                                        Delta_Q         = Delta_Q         ,&
                                        bet             = bet             ,&
                                        nTau            = nTau            ,&
                                        nStages         = nStages         ,&
                                        nI              = nI              ,&
                                        Newtonl2Res     = Newtonl2Res(nI), &
                                        Source_Fac      = Source_Fac      ,&
                                        gam             = gam             ,&
                                        Mach            = Mach            ,&
                                        u               = u               ,&
                                        p               = p               ,&
                                        c               = c               ,&
                                        gm1             = gm1             ,&
                                        rho             = rho             ,&
                                        rho_Static      = rho_Static      ,&
                                        P_Static        = P_Static        ,&
                                        p0In            = p0In            ,&
                                        rho0In          = rho0In          ,&
                                        Q_Exact         = Q_Exact         ,&
                                        Jac_Curv        = Jac_Curv        ,&
                                        Inv_Jac_Curv    = Inv_Jac_Curv    ,&
                                        dZidX           = dZidX           ,&
                                        dZidt           = dZidt           ,&
                                        nFlow           = nFlow           ,&
                                        Area            = Area)             
                        ELSE IF (nL == 2) THEN 
                            CALL TriDi(  iMin           = iMin              ,&
                                         iMax           = iMax              ,&
                                         bMin           = bMin              ,&
                                         bMax           = bMax              ,&
                                         Delta_Q        = Delta_Q           ,&
                                         A_Jac          = A_Jac             ,&
                                         S_Jac          = S_Jac             ,&
                                         Q_n            = Q_n               ,&
                                         Q_np1_k        = Q_np1_k           ,&
                                         Q_np1_kp1      = Q_np1_kp1         ,&
                                         Q_np1          = Q_np1             ,&
                                         Q_nm1          = Q_nm1             ,&
                                         RHS            = RHS               ,&
                                         time           = time              ,&
                                         delta_x        = delta_x           ,&
                                         delta_Zi       = delta_Zi          ,&
                                         delta_tau      = delta_tau         ,&
                                         delta_t        = delta_t           ,&
                                         DS             = DS                ,&
                                         nD             = nD                ,&
                                         Dis            = Dis               ,&
                                         bet            = bet               ,&
                                         nTau           = nTau              ,&
                                         nStages        = nStages           ,&
                                         Scaling_Fac    = Scaling_Fac       ,&
                                         E              = E                 ,&
                                         E_GPT          = E_GPT             ,&
                                         Newtonl2Res    = Newtonl2Res(nI)   ,&
                                         dQdTau         = dQdTau            ,&
                                         Source_Fac     = Source_Fac        ,&
                                         gam            = gam               ,&
                                         Mach           = Mach              ,&
                                         u              = u                 ,&
                                         p              = p                 ,&
                                         c              = c                 ,&
                                         gm1            = gm1               ,&
                                         rho            = rho               ,&
                                         rho_Static     = rho_Static        ,&
                                         P_Static       = P_Static          ,&
                                         p0In           = p0In              ,&
                                         rho0In         = rho0In            ,&
                                         Q_Exact        = Q_Exact           ,&
                                         Jac_Curv       = Jac_Curv          ,&
                                         Inv_Jac_Curv   = Inv_Jac_Curv      ,&
                                         dZidX          = dZidX             ,&
                                         dZidt          = dZidt             ,&
                                         nFlow          = nFlow             ,&
                                         Area           = Area)
                        ELSE
                            CYCLE
                        END IF
                        
                        IF (nI >= 2) THEN
                            NewtonResidual(nI) = ABS(Newtonl2Res(nI) - &
                                            1.0_rDef*Newtonl2Res(nI-1))
                            ELSE
                            NewtonResidual(nI) = Newtonl2Res(nI)
                        END IF

                        !WRITE(0, *) LeftHandSide, DifferStencil,        &
                        !            nT, nI, NewtonResidual(nI),         &
                        !            FINDLOC(dQdTau, MAXVAL(dQdTau)),    &
                        !            MAXVAL(dQdTau)
                         
                        IF (NewtonResidual(nI) > tol2) THEN
                            Q_np1_k = Q_np1_kp1
                        ELSE
                            EXIT
                        END IF
                    END DO

                    Q_n     = Q_np1_kp1
                    
                    l2Fac = 0.0_rDef
                    DO j = bMin, bMax
                        DO i = iMin, iMax 
                            l2Fac = l2Fac+Q_n(i, j)*Q_n(i, j)
                        END DO
                    END DO
                    l2Res(nT)       = SQRT(l2Fac)/REAL((iMax)*3, rDef)
                    
                    IF (nT >= 2) THEN
                        Residual(nT)    = ABS(l2Res(nT-1) - l2Res(nT))
                        ELSE
                        Residual(nT) = l2Res(nT)
                    END IF
                    
                    CALL CreateObject(object = thisCFLObject, & 
                        Residual        = Residual(nT) ,&
                        numTimeSteps    = numTimeSteps     ,&
                        nT              = nT)

                    CALL CFL_Ramping(   object          = thisCFLObject    ,&
                        CFL              = nu_physical  )

                    CALL DestroyObject(object = thisCFLObject)

!
!                    CALL CFL_Ramping(   Residual        = Residual(nT)     ,&
!                                        nu              = nu_physical      ,&
!!
                    !IF (nT == 1) THEN
                    !    CALL CFL_Ramping2(  Residual1   = Residual(nT)  ,&
                    !                        Residual2   = Residual(nT)  ,&
                    !                        nu          = nu_physical)
                    !ELSE
                    !    CALL CFL_Ramping2(  Residual1   = Residual(nT)      ,&
                    !                        Residual2   = Residual(nT-1)    ,&
                    !                        nu          = nu_physical)
                    !END IF

                    IF (MOD(nT, 500) == 501) THEN
                        OPEN(57, FILE = '/Users/zaidhs/Documents/PhD/&
                                    &Third_CAA_C1_P1/TimeStepSol/'&
                                        //filename3, FORM = 'FORMATTED')
                        DO i = iMin,   iMax
                            WRITE(57, *)    x(i)                            ,&
                                            !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                                            Q_n(i, 1)*Jac_Curv(i)           ,&
                                            !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                                            Q_n(i, 2)/Q_n(i, 1)             ,&
                                            !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                                            ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                            - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                            /Q_n(i, 1)))*Jac_Curv(i)        ,&
                                            !!!!!!!!!!!!Mach Number!!!!!!!!!!!
                                            (Q_n(i, 2)/Q_n(i, 1))/(SQRT(gam*(&
                                            ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                            - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                            /Q_n(i, 1))))/(Q_n(i, 1))))     ,& 
                                            !!!!!!!!!!!!Exact Density!!!!!!!!!
                                            Q_Exact(i, 1)*Jac_Curv(i)       ,&
                                            !!!!!!!!!!!!Exact Velocity!!!!!!!!
                                            Q_Exact(i, 2)/Q_Exact(i, 1)     ,&
                                            !!!!!!!!!!!!Exact Pressure!!!!!!!!
                                            ((gam-1.0_rDef)*(Q_Exact(i, 3) &
                                            - 0.5_rDef*Q_Exact(i, 2)*&
                                            Q_Exact(i, 2)/Q_Exact(i, 1)))*&
                                            Jac_Curv(i)                     ,&
                                            !!!!!!!!!!!!Exact Mach Number!!!!!
                                            (Q_Exact(i, 2)/Q_Exact(i, 1))/&
                                            (SQRT(gam*(((gam-1.0_rDef)*&
                                            (Q_Exact(i, 3) - 0.5_rDef*&
                                            Q_Exact(i, 2)*Q_Exact(i, 2)&
                                            /Q_Exact(i, 1))))/(Q_Exact(i, 1)))), &
                                            i, delta_xNew(i), Jac_Curv(i)
                        END DO
                        CLOSE(57)
                    END IF
                    
                    IF (MOD(nT, 1) == 0) THEN
                        WRITE(6, *) 'Mean ', LeftHandSide, ' ',  &
                                    DifferStencil, ' ',  &
                                    Dissipation, nT, nI, Residual(nT), & 
                                    nu_physical, delta_t
                    END IF
        
                    WRITE(56, *)    LeftHandSide, DifferStencil, &
                                    Dissipation, nT, nI, Residual(nT), & 
                                    nu_physical                    

                    IF (Residual(nT) < tol) THEN
                        WRITE(6, *) 'Mean ', LeftHandSide, ' ',  Dissipation, &
                                    ' ', DifferStencil, &
                                    "  Converged After: ", nT, &
                                    "TimeSteps"
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        WRITE(filename7, '(a18)') 'Initialize_CAA.dat'
                        OPEN(68, FILE = '/Users/zaidhs/Documents/PhD/&
                                    &Third_CAA_C1_P1/FinalSol/'&
                                        //filename7, FORM = 'FORMATTED')
                        WRITE(68, *)    p_Exit   ,&
                                        iMax    ,&
                                        iMin    ,&
                                        bMin    ,&
                                        bMax    ,&
                                        gam     ,&
                                        nD      ,&
                                        nL      ,&
                                        DS
                        CLOSE(68)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        WRITE(filename4, '(a17)') 'Mean_Flow_CAA.dat'
                        OPEN(66, FILE = '/Users/zaidhs/Documents/PhD/&
                                    &Third_CAA_C1_P1/FinalSol/'&
                                        //filename4, FORM = 'FORMATTED')
                        DO i = iMin, iMax
                            WRITE(66, *)    Q_n(i, 1)       ,&
                                            Q_n(i, 2)       ,&
                                            Q_n(i, 3)       ,&
                                            Q_Exact(i, 1)   ,&
                                            Q_Exact(i, 2)   ,&
                                            Q_Exact(i, 3)   ,&
                                            x(i)            ,&
                                            Jac_Curv(i)
                        END DO
                        CLOSE(66)
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        EXIT
                    ELSE IF (Residual(nT) > 100) THEN
                        WRITE(6, *) 'Mean ', DifferStencil, "Blew Up: ", &
                                    nT, "TimeSteps, Press Any TO Continue"&
                                    ,Residual(nT)
                        EXIT
                    END IF
                END DO
                CLOSE(56)
            END DO
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO 

    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_Mean(i, j)   = Q_n(i, j)
        END DO
    END DO

     
    END SUBROUTINE MeanValues

END MODULE GetMeanValues
