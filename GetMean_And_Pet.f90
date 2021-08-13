MODULE GetMean_And_Pet
    USE,  INTRINSIC:: ISO_FORTRAN_ENV

    USE GetLUSGS
    USE GetJacobians
    USE GetTriDi
    USE GetAreaDerivative
    USE GetDX
    USE GetExactSol
    USE GetiMaxValue
    USE Precision_Def
    USE FFT1D

    IMPLICIT NONE
    PRIVATE
    
    PUBLIC:: Mean_And_Pet

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE Mean_And_Pet(Q_Mean_And_Pet  ,&
                            Q_Mean          ,&
                            Q_Peturbation   ,&
                            Q_Exact         ,&
                            nFlow           ,&
                            p_Exit          ,&
                            nL              ,&
                            nD              ,&
                            DS)

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: Q_Mean

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(INOUT):: Q_Mean_And_Pet      ,&
                                                        Q_Exact             ,&
                                                        Q_Peturbation

    INTEGER, INTENT(IN)::   nFlow   ,&
                            nL      ,&
                            nD      ,&
                            DS

    REAL(KIND = rDef), INTENT(IN):: p_Exit

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
                        numOfCycles     ,&
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
                        xC              ,&
                        Steps_Per_Cycle

    REAL(KIND = rDef) ::    delta_x         ,&
                            delta_x1        ,&
                            delta_x2        ,&
                            delta_x3        ,&
                            pi              ,&
                            delta_t         ,&
                            delta_Zi        ,&
                            delta_tau       ,&
                            tol             ,&
                            tol2            ,&
                            gam             ,&
                            gm1             ,&
                            omega1          ,&
                            omega2          ,&
                            time            ,&
                            L2Fac           ,&
                            NewtonL2Fac     ,&
                            epsi_A          ,&
                            epsi_S          ,&
                            Scaling_Fac     ,&
                            kDeltaX         ,&
                            sigma_Inv       ,&
                            nu_physical     ,&
                            nu              ,&
                            xx              ,&
                            rho             ,&
                            u               ,&
                            p               ,&
                            Mach            ,&
                            c               ,&
                            eig1_A          ,&
                            eig2_A          ,&
                            eig3_A          ,&
                            eig1_S          ,&
                            eig2_S          ,&
                            eig3_S          ,&
                            T0In            ,&
                            p0In            ,&
                            rho0In          ,&
                            P_Static        ,&
                            rho_Static      ,&
                            FinalTime       ,&
                            Cycle_Time

                    
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
                                                        dAdZi           ,&
                                                        Source_Fac      ,&
                                                        dXdt            ,&
                                                        dZidX           ,&
                                                        dZidt           ,&
                                                        Vel             ,&
                                                        Vel_Mean        ,&
                                                        Vel_Exact       ,&
                                                        Vel_Pet         ,&
                                                        Vel_Exact_Pet   ,&
                                                        Den             ,&
                                                        Den_Mean        ,&
                                                        Den_Exact       ,&
                                                        Den_Pet         ,&
                                                        Den_Exact_Pet   ,&
                                                        Pres            ,&
                                                        Pres_Mean       ,&
                                                        Pres_Exact      ,&
                                                        Pres_Pet        ,&
                                                        Pres_Exact_Pet  ,&
                                                        Pressure_FFT    ,&
                                                        Time_FFT

    REAL(KIND = rDef), DIMENSION(:, :), ALLOCATABLE ::  Q_n                 ,&
                                                        Q_np1               ,&
                                                        Q_nm1               ,&
                                                        dQdTau              ,&
                                                        RHS                 ,&
                                                        E_GPT               ,&
                                                        E                   ,&
                                                        Q_np1_EGPT          ,&
                                                        Q_EGPT              ,&
                                                        Delta_Q_EGPT        ,&
                                                        Q_np1_k             ,&
                                                        Q_np1_kp1           ,&
                                                        Delta_Q_star2       ,&
                                                        Delta_Q_star        ,&
                                                        Delta_Q             ,&
                                                        AA                  ,&
                                                        BB                  ,&
                                                        EE                  ,&
                                                        Dis                 ,&
                                                        Exact               ,&
                                                        Q_Peturb_Out        ,&
                                                        Min_Pressure        ,&
                                                        Min_Pressure_Exact  ,&
                                                        Max_Pressure        ,&
                                                        Max_Pressure_Exact

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
    CHARACTER(LEN = 45):: filename3
    CHARACTER(LEN = 25):: filename4
    CHARACTER(LEN = 38):: filename5
    CHARACTER(LEN = 28):: filename6
    CHARACTER(LEN = 28):: filename7

    pi              = 4.0_rDef*ATAN(1.0_rDef)
    xMin            = -10                 
    xMax            =  10               
    delta_Zi        = 1.0_rDef
    iMin            = 1
    bMin            = 1
    bMax            = 3
    tol             = 1E-12_rDef
    tol2            = 1E-16_rDef
    gam             = 1.40_rDef
    gm1             = 1.0_rDef/(gam-1.0_rDef)
    MaxIteration    = 151
    omega1          = pi/10.0_rDef
    omega2          = pi/50.0_rDef
    time            = 0.0_rDef
    nStages         = 4
    delta_t         = 0.0_rDef
    Cycle_Time      = (2.0_rDef*pi)/(0.6_rDef*pi)
    kDeltaX         = 0.0_rDef  ! Just Needed to Iitialize

    ! Where does this come from and Why is it needed? The Inflow/Outflow 
    ! boundary is unsteady now that I'm solving CAA problem. The
    ! Unsteady boundary claims epsilon*SIN(omega*B + omega*time)) 
    ! or epsilon*COS(omega*V +omega*time)). To find the period of a SIN/COS wave
    ! (or as i'm calling it Cycle_Time) we divide 2pi by omega. Note that 
    ! we dont care about 'B' since that only shifts the
    ! functions so no big deal. Nor is epislon important for the period (cycle)
    ! since that deals with the amplitude of the wave. Omega = 0.6*pi 

    numOfCycles     = 90
    Steps_Per_Cycle = 2**8
    FinalTime       = Cycle_Time*REAL(numOfCycles, rDef)
    numTimeSteps    = INT(Steps_Per_Cycle*INT(numOfCycles)) + 1
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

    ALLOCATE(   Area(iMax)                              ,&
                Area1(iMax)                             ,&
                x1(iMax)                                ,&
                Source_Fac(iMax)                        ,&
                Q_n(iMax, bMax)                         ,&
                Q_Store(numTimeSteps, iMax, bMax)       ,&
                Pressure_FFT(numTimeSteps)              ,&
                TIme_FFT(numTimeSteps)                  ,&
                Q_Peturb_Out(numTimeSteps, bMax)        ,&
                Q_np1(iMax, bMax)                       ,&
                Q_nm1(iMax, bMax)                       ,&
                dQdTau(iMax, bMax)                      ,&
                dAdZi(iMax)                             ,&
                RHS(iMax, bMax)                         ,&
                l2Res(numTimeSteps)                     ,&
                Residual(numTimeSteps)                  ,&
                MaxWaveSpeed(numTimeSteps)              ,&
                Newtonl2Res(MaxIteration)               ,&
                NewtonResidual(MaxIteration)            ,&
                E_GPT(iMax, bMax)                       ,&
                E(iMax, bMax)                           ,&
                Dis(iMax, bMax)                         ,&
                Q_np1_k(iMax, bMax)                     ,&
                Q_np1_kp1(iMax, bMax)                   ,&
                L_LowerDiag(bMax, bMax, iMax)           ,&
                L_Diag(bMax, bMax, iMax)                ,&
                U_UpperDiag(bMax, bMax, iMax)           ,&
                U_Diag(bMax, bMax, iMax)                ,&
                D_Diag(bMax, bMax, iMax)                ,&
                Exact(iMax, bMax)                       ,&
                Delta_Q_star2(iMax, bMax)               ,&
                Delta_Q_star(iMax, bMax)                ,&
                Delta_Q(iMax, bMax)                     ,&
                AA(bMax, bMax)                          ,&
                BB(bMax, bMin)                          ,&
                EE(bMax, bMin)                          ,&
                bet(nStages)                            ,&
                A_Jac(bMax, bMax, iMax)                 ,&
                A_Plus(bMax, bMax, iMax)                ,&
                A_Minus(bMax, bMax, iMax)               ,&
                S_Jac(bMax, bMax, iMax)                 ,&
                S_Plus(bMax, bMax, iMax)                ,&
                S_Minus(bMax, bMax, iMax)               ,&
                Jac_Curv(iMax)                          ,&
                Inv_Jac_Curv(iMax)                      ,&
                Vel(iMax)                               ,&
                Vel_Mean(iMax)                          ,&
                Vel_Exact(iMax)                         ,&
                Vel_Pet(iMax)                           ,&
                Vel_Exact_Pet(iMax)                     ,&
                Den(iMax)                               ,&
                Den_Mean(iMax)                          ,&
                Den_Exact(iMax)                         ,&
                Den_Pet(iMax)                           ,&
                Den_Exact_Pet(iMax)                     ,&
                Pres(iMax)                              ,&
                Pres_Mean(iMax)                         ,&
                Pres_Exact(iMax)                        ,&
                Pres_Pet(iMax)                          ,&
                Pres_Exact_Pet(iMax)                    ,&
                dXdt(iMax)                              ,&
                dZidX(iMax)                             ,&
                dZidt(iMax)                             ,&
                Min_Pressure(numTimeSteps, iMax)        ,&
                Min_Pressure_Exact(numTimeSteps, iMax)  ,&
                Max_Pressure(numTimeSteps, iMax)        ,&
                Max_Pressure_Exact(numTimeSteps, iMax))

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
    
!! RK Coefficients
    bet(1)              = (1.0_rDef)/(4.0_rDef) 
    bet(2)              = (1.0_rDef)/(3.0_rDef) 
    bet(3)              = (1.0_rDef)/(2.0_rDef) 
    bet(4)              =  1.0_rDef 
    
!! Different LHS. 1--> Uses Lower upper Symmetric Gauss Seidel Method by 
!!                      Yoon and Jameson
!!                2--> Uses ADI on the LHS
    IF (nL == 1) THEN
        LeftHandSide = "LUSGS"
    ELSE IF (nL == 2) THEN
        LeftHandSide = "TriDi"
    END IF
!! Different RHS Dissipation.   1--> Second Order Dissipation
!!                              2--> Fourth Order
!!                              3--> Sixth Order
!!                              4--> 8th Order
!!                              5--> Tenth Order
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
    END IF
!! Different RHS Differncing Stencil.   1--> Second Order
!!                                      2--> Fourth Order
!!                                      3--> Sixth Order
!!                                      4--> RDRP
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
    END IF
    
!! Define Area and Domain

    DO i = iMin, iMax
        IF (x(i) >= 0.0_rDef) THEN
            Area(i) = 0.536572_rDef-0.198086_rDef*&
                      EXP(-LOG(2.0_rDef)*((x(i)/0.6_rDef)**2.0_rDef))
        ELSEIF (x(i) < 0.0_rDef) THEN
            Area(i) = 1.0_rDef-0.661514_rDef*&
                      EXP(-LOG(2.0_rDef)*((x(i)/0.6_rDef)**2.0_rDef))
        END IF
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

!WRITE(59, *) nPts
    DO i = iMin, iMax
        Source_Fac(i)  = -Jac_Curv(i)*(dAdZi(i))
    END DO
!! Q_Mean_And_Pet**AT THIS POINT**is Q_mean. We use Q_Mean as an 
!! inital condition. 
    DO j = bMin, bMax
        DO i = iMin, iMax
            Q_n(i, j) = Q_Mean(i, j)
        END DO
    END DO

!! Needed For Outflow Boundary Calculation
    P_Static    = p_Exit*Inv_Jac_Curv(1)
    rho_Static  = Q_n(iMin, 1) 

    WRITE(filename6, '(a5, a1, a8, a1, a4, a9)')&
            LeftHandSide, '_', DifferStencil, '_',           &
            Dissipation, '_Time.dat'
    OPEN(61, FILE = filename6, FORM = 'FORMATTED')


    DO nT = 1, numTimeSteps

        delta_t = Cycle_Time/REAL(Steps_Per_Cycle, rDef)

        WRITE(filename3, '(a5, a1, a8, a1, a4, a1, i8.8, a2, a15)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_', nT, 'nT', '_Mean_Petu.dat' 

        WRITE(filename5, '(a5, a1, a8, a1, a4, a1, i8.8, a2, a8)')&
                LeftHandSide, '_', DifferStencil, '_',           &
                Dissipation, '_', nT, 'nT', 'Petu.dat' 

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
            nu              = 2.5_rDef
            
            IF (MOD(nt, 5) == 1) THEN
                !delta_t         = nu_physical/SQRT(sigma_Inv&
                !                                  *sigma_Inv)
                nu_physical     = delta_t*SQRT(sigma_Inv*sigma_Inv)
                delta_tau       = (nu)/SQRT(sigma_Inv*sigma_Inv &
                                    + (1.0_rDef/delta_t)&
                                    *(1.0_rDef/delta_t))
            END IF

            IF (nL == 1) THEN
                CALL LUSGS( iMin          = iMin            ,&
                            iMax          = iMax            ,&
                            bMin          = bMin            ,&
                            bMax          = bMax            ,&
                            Q_n           = Q_n             ,&
                            Q_np1_k       = Q_np1_k         ,&
                            Q_np1         = Q_np1           ,&
                            Q_nm1         = Q_nm1           ,&
                            Q_np1_kp1     = Q_np1_kp1       ,&
                            RHS           = RHS             ,&
                            time          = time            ,&
                            dQdTau        = dQdTau          ,&
                            delta_x       = delta_x         ,&
                            delta_Zi      = delta_Zi        ,&
                            delta_tau     = delta_tau       ,&
                            delta_t       = delta_t         ,&
                            DS            = DS              ,&
                            nD            = nD              ,&
                            Dis           = Dis             ,&
                            L_LowerDiag   = L_LowerDiag     ,&
                            L_Diag        = L_Diag          ,&
                            U_UpperDiag   = U_UpperDiag     ,&
                            U_Diag        = U_Diag          ,&
                            D_Diag        = D_Diag          ,&
                            A_Plus        = A_Plus          ,&
                            A_Minus       = A_Minus         ,&
                            S_Plus        = S_Plus          ,&
                            S_Minus       = S_Minus         ,&
                            S_Jac         = S_Jac           ,&
                            Delta_Q_star2 = Delta_Q_star2   ,&
                            AA            = AA              ,&
                            BB            = BB              ,&
                            EE            = EE              ,&
                            Delta_Q_star  = Delta_Q_star    ,&
                            Delta_Q       = Delta_Q         ,&
                            bet           = bet             ,&
                            nTau          = nTau            ,&
                            nStages       = nStages         ,&
                            nI            = nI              ,&
                            Newtonl2Res   = Newtonl2Res(nI), &
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

        DO i = iMin, iMax
            DO j = bMin, bMax 
                Q_n(i, j) = Q_np1_kp1(i, j)
            END DO
        END DO

        DO i = iMin, iMax
        !!!!!!!!!!!!!!!!!! Velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Vel(i)              = Q_n(i, 2)/Q_n(i, 1) 
            Vel_Mean(i)         = Q_Mean(i, 2)/Q_Mean(i, 1)
            Vel_Pet(i)          = Vel(i) - Vel_Mean(i)
            Vel_Exact(i)        = Q_Exact(i, 2)/Q_Exact(i, 1)
            Vel_Exact_Pet(i)    = Vel(i) - Vel_Exact(i) 
        !!!!!!!!!!!!!!!!!! Density   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Den(i)              = Q_n(i, 1)*Jac_Curv(i)
            Den_Mean(i)         = Q_Mean(i, 1)*Jac_Curv(i)
            Den_Pet(i)          = Den(i) - Den_Mean(i)
            Den_Exact(i)        = Q_Exact(i, 1)*Jac_Curv(i)
            Den_Exact_Pet(i)    = Den(i) - Den_Exact(i) 
        !!!!!!!!!!!!!!!!!! Pressure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            Pres(i)             = ((gam-1.0_rDef)*(Q_n(i, 3)     &
                                    - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)&
                                    /Q_n(i, 1)))*Jac_Curv(i)
            Pres_Mean(i)        = ((gam-1.0_rDef)*(Q_Mean(i, 3)     &
                                    - 0.5_rDef*Q_Mean(i, 2)*Q_Mean(i, 2)&
                                    /Q_Mean(i, 1)))*Jac_Curv(i)
            Pres_Pet(i)         = Pres(i) - Pres_Mean(i)
            Pres_Exact(i)       = ((gam-1.0_rDef)*(Q_Exact(i, 3)     &
                                    - 0.5_rDef*Q_Exact(i, 2)*Q_Exact(i, 2)&
                                    /Q_Exact(i, 1)))*Jac_Curv(i)
            Pres_Exact_Pet(i)    = Pres(i) - Pres_Exact(i) 
        END DO

        DO i = iMin, iMax
            Max_Pressure(nT, i)         = Pres_Pet(i)
            Min_Pressure(nT, i)         = Pres_Pet(i)
            Max_Pressure_Exact(nT, i)   = Pres_Exact_Pet(i)
            Min_Pressure_Exact(nT, i)   = Pres_Exact_Pet(i)
        END DO
        
        time = time+delta_t
        
        IF (MOD(nT, 1) == 0) THEN
            WRITE(6, *) 'Mean and Petrubation ',                &
                        LeftHandSide, ' ',                      &
                        DifferStencil, ' ',                     &
                        Dissipation, nT, nI,                    & 
                        time, INT(time/Cycle_Time+1.0_rDef),    &
                        nu_physical 
        END IF

        IF (MOD(nT, INT(Steps_Per_Cycle/8)) == 1) THEN
            WRITE(6, *) "Writing Out Solution Here. Check File"
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
                                /Q_Exact(i, 1))))/&
                                (Q_Exact(i, 1)))), &
                                i, delta_xNew(i), Jac_Curv(i)
            END DO
            CLOSE(57)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            OPEN(58, FILE = '/Users/zaidhs/Documents/PhD/&
                        &Third_CAA_C1_P1/TimeStepSol/'&
                            //filename5, FORM = 'FORMATTED')
            DO i = iMin,   iMax
                WRITE(58, *) x(i)                               ,&
                            !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                            Den_Pet(i)                          ,&
                            !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                            Vel_Pet(i)                          ,&
                            !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                            Pres_Pet(i)                         ,&
                            !!!!!!!!!!!!Exact Density!!!!!!!!!
                            Den_Exact_Pet(i)                    ,&
                            !!!!!!!!!!!!Exact Velocity!!!!!!!!
                            Vel_Exact_Pet(i)                    ,&
                            !!!!!!!!!!!!Exact Pressure!!!!!!!!
                            Pres_Exact_Pet(i)                   ,&
                            time, INT(time/Cycle_Time+1.0_rDef), &
                            MAXVAL(Max_Pressure(:, i))          ,&
                            MAXVAL(Max_Pressure_Exact(:, i))    ,&
                            MINVAL(Min_Pressure(:, i))          ,&
                            MINVAL(Min_Pressure_Exact(:, i))    ,&
                            0.5_rDef*(MAXVAL(Max_Pressure(:, i)) &
                            + MINVAL(Min_Pressure(:, i)))
            END DO
            CLOSE(58)
        END IF
        !Solution is written out at every time step 
        WRITE(61, *)    nT                                  ,& 
                        !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                        Den_Pet(iMax)                       ,&
                        !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                        Vel_Pet(iMax)                       ,&
                        !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                        Pres_Pet(iMax)                         ,&
                        !!!!!!!!!!!!Exact Density!!!!!!!!!
                        Den_Exact_Pet(iMax)                    ,&
                        !!!!!!!!!!!!Exact Velocity!!!!!!!!
                        Vel_Exact_Pet(iMax)                    ,&
                        !!!!!!!!!!!!Exact Pressure!!!!!!!!
                        Pres_Exact_Pet(iMax)                   ,&
                        time                                ,&
                        INT(time/Cycle_Time+1.0_rDef)

        Pressure_FFT(nT)    = Pres_Pet(iMax) 
        Time_FFT(nT)        = time

        IF (ABS(FinalTime-time) <= tol) THEN
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
                                /Q_Exact(i, 1))))/&
                                (Q_Exact(i, 1)))), &
                                i, delta_xNew(i), Jac_Curv(i), &
                                MAXVAL(Max_Pressure(:, i))          ,&
                                MAXVAL(Max_Pressure_Exact(:, i))    ,&
                                MINVAL(Min_Pressure(:, i))          ,&
                                MINVAL(Min_Pressure_Exact(:, i))    ,&
                                0.5_rDef*(MAXVAL(Max_Pressure(:, i)) &
                                + MINVAL(Min_Pressure(:, i)))
            END DO
            CLOSE(57)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            OPEN(58, FILE = '/Users/zaidhs/Documents/PhD/&
                        &Third_CAA_C1_P1/TimeStepSol/'&
                            //filename5, FORM = 'FORMATTED')
            DO i = iMin,   iMax
                WRITE(58, *) x(i)                               ,&
                            !!!!!!!!!!!!Density!!!!!!!!!!!!!!!
                            Den_Pet(i)                          ,&
                            !!!!!!!!!!!!Velocity!!!!!!!!!!!!!!
                            Vel_Pet(i)                          ,&
                            !!!!!!!!!!!!Pressure!!!!!!!!!!!!!!
                            Pres_Pet(i)                         ,&
                            !!!!!!!!!!!!Exact Density!!!!!!!!!
                            Den_Exact_Pet(i)                    ,&
                            !!!!!!!!!!!!Exact Velocity!!!!!!!!
                            Vel_Exact_Pet(i)                    ,&
                            !!!!!!!!!!!!Exact Pressure!!!!!!!!
                            Pres_Exact_Pet(i)                   ,&
                            time, INT(time/Cycle_Time+1.0_rDef), &
                            MAXVAL(Max_Pressure(:, i))          ,&
                            MAXVAL(Max_Pressure_Exact(:, i))    ,&
                            MINVAL(Min_Pressure(:, i))          ,&
                            MINVAL(Min_Pressure_Exact(:, i))    ,&
                            0.5_rDef*(MAXVAL(Max_Pressure(:, i)) &
                            + MINVAL(Min_Pressure(:, i)))
            END DO
            CLOSE(58)
        END IF
    END DO
    CLOSE(61)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL Real_FFT(  data    = Pressure_FFT(1:Steps_Per_Cycle)   ,&
                    n       = Steps_Per_Cycle                   ,&
                    iSign   = 1_i_def)

    WRITE(filename7, '(a5, a1, a8, a1, a4, a9)')&
            LeftHandSide, '_', DifferStencil, '_',           &
            Dissipation, '_FFFT.dat'
    OPEN(72, FILE = filename7, FORM = 'FORMATTED')
    
    DO nT = 1, Steps_Per_Cycle
        WRITE(72, *) nT, Time_FFT(nT), Pressure_FFT(nT)
    END DO
    CLOSE(72)
     
    END SUBROUTINE Mean_And_Pet

END MODULE GetMean_And_Pet
