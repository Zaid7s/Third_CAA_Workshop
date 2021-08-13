MODULE GetJacobians
!! Working TestCase

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: Jacobians

    INTEGER, PARAMETER:: rDef = REAL64

    CONTAINS

    SUBROUTINE Jacobians(  iMin         ,&
                           iMax         ,&
                           bMin         ,&
                           bMax         ,&
                           A_Jac        ,&
                           u            ,&
                           p            ,&
                           c            ,&
                           rho          ,&
                           Mach         ,&
                           Q_n          ,&
                           gm1          ,&
                           gam          ,&
                           eig1_A       ,&
                           eig2_A       ,&
                           eig3_A       ,&
                           epsi_A       ,&
                           Source_Fac   ,&
                           A_Plus       ,&
                           A_Minus      ,&
                           eig1_S       ,&
                           eig2_S       ,&
                           eig3_S       ,&
                           epsi_S       ,&
                           S_Plus       ,&
                           S_Minus      ,& 
                           S_Jac        ,&
                           Jac_Curv     ,&
                           Inv_Jac_Curv, &
                           dZidX        ,&
                           dZidt        ,&
                           Area)
    
    INTEGER, INTENT (IN):: iMin        ,&
                            iMax        ,&
                            bMin        ,&
                            bMax

    REAL(KIND = rDef), INTENT(INOUT)   :: eig1_A  ,&
                                        eig2_A  ,&
                                        eig3_A  ,&
                                        epsi_A  ,&
                                        eig1_S  ,&
                                        eig2_S  ,&
                                        eig3_S  ,&
                                        epsi_S  
        
    REAL(KIND = rDef), INTENT(IN)    :: gm1     ,&
                                        gam

    REAL(KIND = rDef), INTENT(INOUT):: u       ,&
                                        p       ,&
                                        c       ,&
                                        Mach    ,&
                                        rho
                                        
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) ::  Jac_Curv        ,&
                                                    Inv_Jac_Curv    ,&
                                                    dZidX           ,&
                                                    dZidt           ,&
                                                    Source_Fac      ,&
                                                    Area

    REAL(KIND = rDef), DIMENSION(:, :), INTENT(IN):: Q_n

    REAL(KIND = rDef), DIMENSION(iMax) ::   Store_EigA      ,&
                                            Store_EigS      ,&
                                            AreaFac

    REAL(KIND = rDef), DIMENSION(:, :, :), INTENT(INOUT):: A_Jac       ,&
                                                            A_Plus      ,&
                                                            A_Minus     ,&
                                                            S_Jac       ,&
                                                            S_Plus      ,&
                                                            S_Minus
    INTEGER:: i, j
    
    DO i = iMin, iMax
        AreaFac(i) = Source_Fac(i)/Area(i)
    END DO

    DO i = iMin, iMax
        A_Jac(1, 1, i) = 0.0_rDef
        A_Jac(1, 2, i) = dZidX(i)*1.0_rDef
        A_Jac(1, 3, i) = 0.0_rDef

        A_Jac(2, 1, i) = dZidX(i)*(Q_n(i, 2)*Q_n(i, 2)*(gam-3.0_rDef))/&
                         (2.0_rDef*Q_n(i, 1)*Q_n(i, 1))
        A_Jac(2, 2, i) = -dZidX(i)*(Q_n(i, 2)*(gam-3.0_rDef))/Q_n(i, 1)
        A_Jac(2, 3, i) = dZidX(i)*(gam-1.0_rDef)
        
        A_Jac(3, 1, i) = -dZidX(i)*((Q_n(i, 2)*(Q_n(i, 2)*Q_n(i, 2)   -   &
                                 Q_n(i, 2)*Q_n(i, 2)*gam +   &
                                 Q_n(i, 1)*Q_n(i, 3)*gam))/  &
                                (Q_n(i, 1)*Q_n(i, 1)*Q_n(i, 1)))
        A_Jac(3, 2, i) = dZidX(i)*(3.0_rDef*Q_n(i, 2)*Q_n(i, 2)      -   &
                       3.0_rDef*Q_n(i, 2)*Q_n(i, 2)*gam  +   &
                       2.0_rDef*Q_n(i, 1)*Q_n(i, 3)*gam)/    &
                       (2.0_rDef*Q_n(i, 1)*Q_n(i, 1))
        A_Jac(3, 3, i) = dZidX(i)*(Q_n(i, 2)*gam)/(Q_n(i, 1))
    END DO
    
    DO i = iMin, iMax
        rho     = Q_n(i, 1)*Jac_Curv(i)
        u       = Q_n(i, 2)/Q_n(i, 1)
        p       = (gam-1.0_rDef)*(Q_n(i, 3) - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)/&
                        Q_n(i, 1))*Jac_Curv(i)
        c       = SQRT(gam*p/rho)
        Mach    = u/c

        eig1_A = u
        eig2_A = u+c
        eig3_A = u-c
        
        Store_EigA(i) = ABS(dZidt(i) + u*dZidX(i)) + c*SQRT(dZidX(i)*dZidX(i))  

    END DO

    epsi_A = MAXVAL(ABS(Store_EigA(:)))
    DO i = iMin, iMax
        A_Plus(1, 1, i)     = 0.5_rDef*(A_Jac(1, 1, i) + epsi_A ) 
        A_Plus(1, 2, i)     = 0.5_rDef*(A_Jac(1, 2, i)          )  
        A_Plus(1, 3, i)     = 0.5_rDef*(A_Jac(1, 3, i)          )  
        A_Plus(2, 1, i)     = 0.5_rDef*(A_Jac(2, 1, i)          )  
        A_Plus(2, 2, i)     = 0.5_rDef*(A_Jac(2, 2, i) + epsi_A )  
        A_Plus(2, 3, i)     = 0.5_rDef*(A_Jac(2, 3, i)          ) 
        A_Plus(3, 1, i)     = 0.5_rDef*(A_Jac(3, 1, i)          ) 
        A_Plus(3, 2, i)     = 0.5_rDef*(A_Jac(3, 2, i)          ) 
        A_Plus(3, 3, i)     = 0.5_rDef*(A_Jac(3, 3, i) + epsi_A )

        A_Minus(1, 1, i)    = 0.5_rDef*(A_Jac(1, 1, i) - epsi_A )
        A_Minus(1, 2, i)    = 0.5_rDef*(A_Jac(1, 2, i)          )  
        A_Minus(1, 3, i)    = 0.5_rDef*(A_Jac(1, 3, i)          )  
        A_Minus(2, 1, i)    = 0.5_rDef*(A_Jac(2, 1, i)          )  
        A_Minus(2, 2, i)    = 0.5_rDef*(A_Jac(2, 2, i) - epsi_A )  
        A_Minus(2, 3, i)    = 0.5_rDef*(A_Jac(2, 3, i)          ) 
        A_Minus(3, 1, i)    = 0.5_rDef*(A_Jac(3, 1, i)          )      
        A_Minus(3, 2, i)    = 0.5_rDef*(A_Jac(3, 2, i)          ) 
        A_Minus(3, 3, i)    = 0.5_rDef*(A_Jac(3, 3, i) - epsi_A ) 
    END DO

    DO i = iMin, iMax
        S_Jac(1, 1, i) = 0.0_rDef
        S_Jac(1, 2, i) = AreaFac(i)*1.0_rDef
        S_Jac(1, 3, i) = 0.0_rDef
        
        S_Jac(2, 1, i) = -AreaFac(i)*(Q_n(i, 2)*Q_n(i, 2))/&
                          (Q_n(i, 1)*Q_n(i, 1))
        S_Jac(2, 2, i) = AreaFac(i)*2.0_rDef*Q_n(i, 2)/Q_n(i, 1)
        S_Jac(2, 3, i) = 0.0_rDef               

        S_Jac(3, 1, i) = -AreaFac(i)*1.0_rDef&
                          *(Q_n(i, 2)*(Q_n(i, 2)*Q_n(i, 2)   -   &
                            Q_n(i, 2)*Q_n(i, 2)*gam          +   &
                            Q_n(i, 1)*Q_n(i, 3)*gam))&
                            /(Q_n(i, 1)*Q_n(i, 1)*Q_n(i, 1))
        S_Jac(3, 2, i) =  AreaFac(i)*(3.0_rDef*Q_n(i, 2)*Q_n(i, 2) -     &
                           3.0_rDef*Q_n(i, 2)*Q_n(i, 2)*gam +               &
                           2.0_rDef*Q_n(i, 1)*Q_n(i, 3)*gam)/&
                           (2.0_rDef*Q_n(i, 1)*Q_n(i, 1))
        S_Jac(3, 3, i) = AreaFac(i)*Q_n(i, 2)*gam/Q_n(i, 1)
    END DO

    DO i = iMin, iMax
        rho     = Q_n(i, 1)*Jac_Curv(i)
        u       = Q_n(i, 2)/Q_n(i, 1)
        p       = (gam-1.0_rDef)*(Q_n(i, 3) - 0.5_rDef*Q_n(i, 2)*Q_n(i, 2)/&
                        Q_n(i, 1))*Jac_Curv(i)
        c       = SQRT(gam*p/rho)
        Mach    = u/c

        eig1_S = u
        eig2_S = u
        eig3_S = u*gam
    
        Store_EigS(i) = gam*ABS(dZidt(i) + u*dZidX(i))
    END DO

    epsi_S = MAXVAL(ABS(Store_EigS(:)))

    DO i = iMin, iMax
        S_Plus(1, 1, i)     = 0.5_rDef*(S_Jac(1, 1, i) + epsi_S ) 
        S_Plus(1, 2, i)     = 0.5_rDef*(S_Jac(1, 2, i)          )  
        S_Plus(1, 3, i)     = 0.5_rDef*(S_Jac(1, 3, i)          )  
        S_Plus(2, 1, i)     = 0.5_rDef*(S_Jac(2, 1, i)          )  
        S_Plus(2, 2, i)     = 0.5_rDef*(S_Jac(2, 2, i) + epsi_S )  
        S_Plus(2, 3, i)     = 0.5_rDef*(S_Jac(2, 3, i)          ) 
        S_Plus(3, 1, i)     = 0.5_rDef*(S_Jac(3, 1, i)          ) 
        S_Plus(3, 2, i)     = 0.5_rDef*(S_Jac(3, 2, i)          ) 
        S_Plus(3, 3, i)     = 0.5_rDef*(S_Jac(3, 3, i) + epsi_S )

        S_Minus(1, 1, i)    = 0.5_rDef*(S_Jac(1, 1, i) - epsi_S )
        S_Minus(1, 2, i)    = 0.5_rDef*(S_Jac(1, 2, i)          )  
        S_Minus(1, 3, i)    = 0.5_rDef*(S_Jac(1, 3, i)          )  
        S_Minus(2, 1, i)    = 0.5_rDef*(S_Jac(2, 1, i)          )  
        S_Minus(2, 2, i)    = 0.5_rDef*(S_Jac(2, 2, i) - epsi_S )  
        S_Minus(2, 3, i)    = 0.5_rDef*(S_Jac(2, 3, i)          ) 
        S_Minus(3, 1, i)    = 0.5_rDef*(S_Jac(3, 1, i)          )      
        S_Minus(3, 2, i)    = 0.5_rDef*(S_Jac(3, 2, i)          ) 
        S_Minus(3, 3, i)    = 0.5_rDef*(S_Jac(3, 3, i) - epsi_S ) 
    END DO

    END SUBROUTINE Jacobians

END MODULE GetJacobians
