MODULE GetExactSol

    USE,  INTRINSIC:: ISO_FORTRAN_ENV
    USE Compressible1DFunctions
    
    IMPLICIT NONE
    PRIVATE
    PUBLIC:: ExactSol

    INTEGER, PARAMETER:: rDef = SELECTED_REAL_KIND(12)

    CONTAINS

    SUBROUTINE ExactSol(  iMin      ,&
                          iMax      ,&
                          p0In      ,&
                          t0In      ,&
                          p_Exit    ,&
                          tol       ,&
                          ii        ,&
                          x         ,&
                          A         ,&
                          p         ,&
                          rho       ,&
                          u)   
    
    INTEGER, INTENT(INOUT) ::   iMin    ,&
                                iMax    ,&
                                ii

    REAL(KIND = rDef), INTENT(INOUT)::  p0In    ,&
                                        t0In    ,&
                                        p_Exit  ,&
                                        tol

    REAL(KIND = rDef), DIMENSION(:), INTENT(INOUT) ::   x       ,&
                                                        A       ,&
                                                        p       ,&
                                                        rho     ,&
                                                        u

    REAL(KIND = rDef):: mIn, p0Out, t0Out, mOut, &
                       aStarIn, aStarOut, aShock, aStar

    REAL(KIND = rDef), PARAMETER:: gam = 1.4_rDef

    REAL(KIND = rDef):: aMin, m1Lo, m1Hi, m1, m2, p0Test, err, t
    
    INTEGER:: aMinLoc, shockLoc

    INTEGER:: i

    REAL(KIND = rDef), DIMENSION(:), ALLOCATABLE:: M,  p0, t0

    LOGICAL:: shockIsPresent, supersonicFlow

    ALLOCATE(M(iMax), &
             p0(iMax), t0(iMax))

    OPEN(16, FILE = 'mach.dat',FORM = 'FORMATTED')

    aMin = MINVAL(A, 1)
    aMinLoc = MINLOC(A, 1)

    WRITE(6, *) 'aMin = ',aMin, aMinLoc

    i = iMax
    CALL GetMFromP0(p0             = p0In,    &
                    p              = p_Exit,    &
                    supersonicFlow = .FALSE., &
                    M              = mOut,    &
                    tol            = tol)

    CALL GetAStar(M     = mOut,      &
                  A     = A(iMax), &
                  AStar = AStarOut)

    WRITE(6, *) 'aStarOut = ',aStarOut

    IF (aStarOut > aMin) THEN
     shockIsPresent = .TRUE.
     aStarIn = aMin
    ELSE
     shockIsPresent = .FALSE.
     aStarIn = aStarOut
    END IF

    WRITE(6, *) 'Is there a shock present? ',shockIsPresent

!   at this point, I can start from the inflow and move to the
!    sonic point (if present)

    supersonicFlow = .FALSE.

    IF (shockIsPresent) THEN

     DO i = 1, aMinLoc
      CALL GetMFromAStar(A              = A(i),           &
                         AStar          = AStarIn,        &
                         supersonicFlow = supersonicFlow, &
                         M              = M(i),           &
                         tol            = tol)      
      p0(i) = p0In
      t0(i) = t0In
      WRITE(6, *) i, x(i), A(i), M(i)
     END DO

!   get the shock data, if needed

     m1Lo = 1.0_rDef
     CALL GetMFromAStar(A              = A(iMax), &
                        AStar          = AStarIn,   &
                        supersonicFlow = .TRUE.,    &
                        M              = m1Hi,      &
                        tol            = tol)

 15 CONTINUE
      m1 = 0.5_rDef*(m1Lo+m1Hi)
      CALL NormalShockRelations(M1  = m1,    &
                                p01 = p0In,  &
                                t01 = t0In,  &
                                M2  = m2,    &
                                p02 = p0Out, &
                                t02 = t0Out)
             
      CALL GetA(M     = m1,      &
                AStar = AStarIn, &
                A     = AShock)

      CALL GetAStar(M     = m2,      &
                    A     = AShock,  &
                    AStar = AStarOut)

      CALL GetMFromAStar(A              = A(iMax), &
                         AStar          = AStarOut,  &
                         supersonicFlow = .FALSE.,   &
                         M              = MOut,      &
                         tol            = tol)      

      CALL GetP0(M  = MOut, &
                 P  = p_Exit, &
                 P0 = p0Test)

      err = p0Out-p0Test

      IF (ABS(err) <= tol) THEN
       CONTINUE
      ELSE
       IF (err > 0.0_rDef) THEN
        m1Lo = m1
       ELSE
        m1Hi = m1
       END IF
       GO TO 15
      END IF

!   shock located.

      WRITE(6, *) 'aShock = ',aShock, m1, m2

      supersonicFlow = .TRUE.
      aStar = aStarIn
      p0    = p0In
      t0    = t0In

      DO i = aMinLoc+1, iMax
       IF (supersonicFlow) THEN
        IF (A(i) > AShock) THEN  ! switch to subsonic
         supersonicFlow = .FALSE.
         AStar = AStarOut
         p0(i) = p0Out
         t0(i) = t0Out
        ELSE
         p0(i) = p0In
         t0(i) = t0In
        END IF
       ELSE
        p0(i) = p0Out
        t0(i) = t0Out
       END IF

       CALL GetMFromAStar(A              = A(i),           &
                          AStar          = AStar,          &
                          supersonicFlow = supersonicFlow, &
                          M              = M(i),           &
                          tol            = tol)      

       WRITE(6, *) i, x(i), A(i), M(i)

      END DO
    ELSE
     DO i = 1, iMax
      CALL GetMFromAStar(A              = A(i),    &
                         AStar          = AStarIn, &
                         supersonicFlow = .FALSE., &
                         M              = M(i),    &
                         tol            = tol)      
      WRITE(6, *) 'Exact Solution', i, x(i), A(i), M(i)
      p0(i) = p0In
      t0(i) = t0In
     END DO
    END IF

!   and write out solution

    OPEN(19, FILE = 'exact.dat',FORM = 'formatted')

    DO i = 1, iMax

     CALL GetP(M  = M(i),  &
               P0 = p0(i), &
               p  = p(i))   

     CALL GetT(M  = M(i),  &
               T0 = T0(i), &
               T  = t)   

     rho(i) = gam*p(i)/t
     u(i)   = M(i)*SQRT(gam*p(i)/rho(i))

     WRITE(16, *) x(i), m(i)
     WRITE(19, *) x(i), p(i), rho(i), u(i)

    END DO

    CLOSE(16)
    CLOSE(19)
    END SUBROUTINE ExactSol

END MODULE GetExactSol
