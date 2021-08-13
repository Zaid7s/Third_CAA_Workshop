MODULE GetDeltaA

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: DeltaA

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS
         

    SUBROUTINE DeltaA(  iMin            ,&
                        iMax            ,&
                        bMin            ,&
                        bMax            ,&
                        delta_tau       ,&
                        delta_t         ,&
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
                        Delta_A_1       ,&
                        Delta_A_2       ,&
                        rho             ,&
                        u               ,&
                        p)

    INTEGER, INTENT(IN):: iMin     ,&
                           iMax     ,&
                           bMax     ,&
                           bMin

    REAL(KIND = rDef), INTENT(IN) ::    gam         ,&
                                        Mean_Mach   ,&
                                        Mean_u      ,&
                                        Mean_p      ,&
                                        Mean_c      ,&
                                        gm1         ,&
                                        Mean_rho    ,&
                                        rho_Static  ,&
                                        P_Static    ,&
                                        p0In        ,&
                                        rho0In      ,&
                                        rho         ,&
                                        u           ,&
                                        p
    
    REAL(KIND = rDef), INTENT(IN):: delta_tau  ,&
                                     delta_t
                                        
    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN) ::   Q_np1_k
    
    REAL(KIND = rDef), INTENT(INOUT) ::   Delta_A_1     ,&
                                        Delta_A_2

    REAL(KIND = rDef) ::    dp0_dp          ,&
                            dp0_du          ,&
                            dp0_drho        ,&
                            dp_dA1          ,&
                            du_dA1          ,&
                            drho_dA1        ,&
                            dp_dA2          ,&
                            du_dA2          ,&
                            drho0_drho      ,&
                            drho0_du        ,&
                            drho0_dp        ,&
                            drho_dA2        ,&
                            p0_current      ,&
                            rho0_current    ,&
                            c               ,&
                            Mach            ,&
                            Delta_p0        ,&
                            Delta_rho0      ,&
                            FacGam1         ,&
                            FacGam2         ,&
                            FacGam3         ,&
                            FacA            ,&
                            FacB            ,&
                            FacC            ,&
                            FacD

    INTEGER:: i, j

    FacGam1  = (gam)/(gam-1.0_rDef)
    FacGam2  = (1.0_rDef)/(gam-1.0_rDef)
    FacGam3  = (2.0_rDef-gam)/(gam-1.0_rDef)

    c       = SQRT(gam*p/rho)
    Mach    = u/c

    p0_current      = p*(1.0_rDef+0.5_rDef*(gam-1.0_rDef)&
                                *Mach*Mach)**FacGam1
    rho0_current    = rho*(1.0_rDef+0.5_rDef*(gam-1.0_rDef)&
                                *Mach*Mach)**FacGam2
    
    Delta_p0    = (p0In     - p0_current    )/delta_t
    Delta_rho0  = (rho0In   - rho0_current  )/delta_t

    dp0_dp      = ((2.0_rDef*gam*p - rho*u*u)*((rho*u*u/(2.0_rDef*p)) &
                    - (rho*u*u/(2.0_rDef*gam*p)) + 1.0_rDef)**(FacGam1))&
                    /(2.0_rDef*gam*p - rho*u*u+gam*rho*u*u)
    dp0_du      = ((rho*u)*((0.5_rDef*rho*u*u*(gam-1.0_rDef))/(gam*p) &
                    + 1.0_rDef)**(FacGam1-1.0_rDef))
    dp0_drho    = ((0.5_rDef*u*u)*((0.5_rDef*rho*u*u*(gam-1.0_rDef))/(gam*p) &
                    + 1.0_rDef)**(FacGam1-1.0_rDef))

    drho0_dp    = -((0.5_rDef*rho*rho*u*u)*(((0.5_rDef*rho*u*u*(gam-1.0_rDef))/(gam*p) &
                    + 1.0_rDef)**(FacGam2-1.0_rDef)))/(gam*p)
    drho0_du    = ((rho*rho*u)*(((0.5_rDef*rho*u*u*(gam-1.0_rDef))/(gam*p) &
                    + 1.0_rDef)**(FacGam2-1.0_rDef)))/(gam*p)
    drho0_drho  = ((gam*((0.5_rDef*rho*u*u*(gam-1.0_rDef)/(gam*p)) & 
                    + 1.0_rDef)**(FacGam2))*(rho*u*u+2.0_rDef*p))/&
                    (2.0_rDef*gam*p - rho*u*u+gam*rho*u*u)
    
    dp_dA1      = 0.0_rDef
    du_dA1      = 0.0_rDef
    drho_dA1    = (-1.0_rDef)/(Mean_c*Mean_c)

    dp_dA2      = (0.5_rDef) 
    du_dA2      = (0.5_rDef)/(Mean_rho*Mean_c)
    drho_dA2    = (0.5_rDef)/(Mean_c*Mean_c)

    FacA        = (dp0_dp*dp_dA1+dp0_du*du_dA1+dp0_drho*drho_dA1)
    FacB        = (dp0_dp*dp_dA2+dp0_du*du_dA2+dp0_drho*drho_dA2)

    FacC        = (drho0_drho*drho_dA1+drho0_du*du_dA1+drho0_dp*dp_dA1)
    FacD        = (drho0_drho*drho_dA2+drho0_du*du_dA2+drho0_dp*dp_dA2)

    Delta_A_1   = (FacB*Delta_rho0-FacD*Delta_p0)/&
                    (FacB*FacC-FacA*FacD)
    Delta_A_2   = -1.0_rDef*(FacA*Delta_rho0-FacC*Delta_p0)/&
                            (FacC*FacB-FacA*FacD)

    END SUBROUTINE DeltaA

END MODULE GetDeltaA
