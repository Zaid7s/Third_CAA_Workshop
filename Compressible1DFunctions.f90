MODULE Compressible1DFunctions
  USE, INTRINSIC:: ISO_FORTRAN_ENV
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: GetMFromP0,    &
            GetMFromAStar, &
            GetP0,         &
            GetAStar,      &
            GetA,          &
            GetP,          &
            GetT,          &
            NormalShockRelations

  INTEGER, PARAMETER:: rDef = REAL64, &
                        r2Def = REAL128
  REAL(KIND = rDef), PARAMETER:: gam = 1.4_rDef

INTERFACE GetMFromP0
  MODULE PROCEDURE GetMFromP0v1
END INTERFACE

INTERFACE GetMFromAStar
  MODULE PROCEDURE GetMFromAStarv1
END INTERFACE

INTERFACE GetAStar
  MODULE PROCEDURE GetAStarFromM
  MODULE PROCEDURE GetAStarFromMQP
END INTERFACE

INTERFACE GetA
  MODULE PROCEDURE GetAFromM
END INTERFACE

INTERFACE GetP
  MODULE PROCEDURE GetPFromM
END INTERFACE

INTERFACE GetT
  MODULE PROCEDURE GetTFromM
END INTERFACE

INTERFACE GetP0
  MODULE PROCEDURE GetP0FromM
  MODULE PROCEDURE GetP0FromMQP
END INTERFACE

INTERFACE NormalShockRelations
  MODULE PROCEDURE NormalShockRelationsv1
END INTERFACE

CONTAINS

SUBROUTINE GetMFromP0v1(p0, p, supersonicFlow, M, tol)
  REAL(KIND = rDef), INTENT(IN):: p0, p, tol
  LOGICAL, INTENT(IN):: supersonicFlow
  REAL(KIND = rDef), INTENT(OUT):: M

! local variables

  REAL(KIND = r2Def):: p0DP, pDP, mDP
  REAL(KIND = r2Def):: mHi, mLo, err, p0test, fac, pow, gm1, gm1o2

  pDP = REAL(p, r2Def)
  p0DP = REAL(p0, r2Def)

  IF (supersonicFlow) THEN
   mHi = 10.0_r2Def
   mLo = 1.0_r2Def
  ELSE
   mHi = 1.0_r2Def
   mLo = 0.0_r2Def
  END IF

 10 CONTINUE
   mDP = 0.5_r2Def*(mHi+mLo)

   CALL GetP0(M  = mDP,    &
              P  = pDP,    &
              P0 = P0Test)

   err = p0DP-p0Test

   IF (ABS(err) <= tol) THEN 
    m = REAL(mDP, rDef)
    RETURN
   ELSE
    IF (err > 0.0_r2Def) THEN
     mLo = mDP
    ELSE
     mHi = mDP
    END IF
   END IF
   GO TO 10
  RETURN 

END SUBROUTINE GetMFromP0v1

SUBROUTINE GetMFromAStarv1(A, AStar, supersonicFlow, M, tol)
  REAL(KIND = rDef), INTENT(IN):: A, AStar, tol
  LOGICAL, INTENT(IN):: supersonicFlow
  REAL(KIND = rDef), INTENT(OUT):: M

! local variables

  REAL(KIND = r2Def):: ADP, AStarDP, mDP
  REAL(KIND = r2Def):: mHi, mLo, err, aSTest, fac, pow, gm1, gm1o2

  aDP = REAL(A, r2Def)
  aStarDP = REAL(AStar, r2Def)
  IF (supersonicFlow) THEN
   mHi = 100.0_r2Def
   mLo = 1.0_r2Def
  ELSE
   mHi = 1.0_r2Def
   mLo = 0.0_r2Def
  END IF

 10 CONTINUE
   mDP = 0.5_r2Def*(mHi+mLo)

   CALL GetAStar(M     = mDP,    &
                 A     = ADP,    &
                 AStar = aSTest)

   err = aStarDP-aSTest

   IF (ABS(err) <= tol) THEN 
    m = REAL(mDP, rDef)
    RETURN
   ELSE
    IF (supersonicFlow) THEN
     IF (err > 0.0_rDef) THEN
      mHi = mDP
     ELSE
      mLo = mDP
     END IF
    ELSE
     IF (err > 0.0_rDef) THEN
      mLo = mDP
     ELSE
      mHi = mDP
     END IF
    END IF
   END IF
   GO TO 10
  RETURN 

END SUBROUTINE GetMFromAStarv1

SUBROUTINE GetP0FromM(M, P, P0)
  REAL(KIND = rDef), INTENT(IN):: M, P
  REAL(KIND = rDef), INTENT(OUT):: P0

! local variables

  REAL(KIND = rDef):: fac, pow, gm1, gp1, gm1o2, aOverAStar2

  gm1 = gam-1.0_rDef
  gm1o2 = 0.5_rDef*gm1

  pow = gam/gm1

  fac = (1.0_rDef+gm1o2*m*m)
  P0 = P*(fac**pow)
  
  RETURN 
END SUBROUTINE GetP0FromM

SUBROUTINE GetP0FromMQP(M, P, P0)
  REAL(KIND = r2Def), INTENT(IN):: M, P
  REAL(KIND = r2Def), INTENT(OUT):: P0

! local variables

  REAL(KIND = r2Def):: fac, pow, gm1, gp1, gm1o2, aOverAStar2

  gm1 = gam-1.0_r2Def
  gm1o2 = 0.5_r2Def*gm1

  pow = gam/gm1

  fac = (1.0_r2Def+gm1o2*m*m)
  P0 = P*(fac**pow)
  
  RETURN 
END SUBROUTINE GetP0FromMQP

SUBROUTINE GetPFromM(M, P0, P)
  REAL(KIND = rDef), INTENT(IN):: M, P0
  REAL(KIND = rDef), INTENT(OUT):: P

! local variables

  REAL(KIND = rDef):: fac, pow, gm1, gp1, gm1o2, aOverAStar2

  gm1 = gam-1.0_rDef
  gm1o2 = 0.5_rDef*gm1

  pow = gam/gm1

  fac = (1.0_rDef+gm1o2*m*m)
  P = P0/(fac**pow)
  
  RETURN 
END SUBROUTINE GetPFromM

SUBROUTINE GetTFromM(M, T0, T)
  REAL(KIND = rDef), INTENT(IN):: M, T0
  REAL(KIND = rDef), INTENT(OUT):: T

! local variables

  REAL(KIND = rDef):: fac, gm1, gp1, gm1o2, aOverAStar2

  gm1 = gam-1.0_rDef
  gm1o2 = 0.5_rDef*gm1

  fac = (1.0_rDef+gm1o2*m*m)
  T = T0/fac
  
  RETURN 
END SUBROUTINE GetTFromM

SUBROUTINE GetAStarFromM(M, A, AStar)
  REAL(KIND = rDef), INTENT(IN):: M, A
  REAL(KIND = rDef), INTENT(OUT):: AStar

! local variables

  REAL(KIND = rDef):: fac, pow, gm1, gp1, gm1o2, aOverAStar2

  gp1 = gam+1.0_rDef
  gm1 = gam-1.0_rDef
  gm1o2 = 0.5_rDef*gm1

  pow = gp1/gm1

  fac = (2.0_rDef/gp1)*(1.0_rDef+gm1o2*m*m)
  aOverAStar2 = (1.0_rDef/(M*M)) &
                *(fac**pow) 
  
  aStar = A/SQRT(aOverAStar2)

  RETURN 
END SUBROUTINE GetAStarFromM

SUBROUTINE GetAStarFromMQP(M, A, AStar)
  REAL(KIND = r2Def), INTENT(IN):: M, A
  REAL(KIND = r2Def), INTENT(OUT):: AStar

! local variables

  REAL(KIND = r2Def):: fac, pow, gm1, gp1, gm1o2, aOverAStar2

  gp1 = gam+1.0_r2Def
  gm1 = gam-1.0_r2Def
  gm1o2 = 0.5_r2Def*gm1

  pow = gp1/gm1

  fac = (2.0_r2Def/gp1)*(1.0_r2Def+gm1o2*m*m)
  aOverAStar2 = (1.0_r2Def/(M*M)) &
                *(fac**pow) 
  
  aStar = A/SQRT(aOverAStar2)

  RETURN 
END SUBROUTINE GetAStarFromMQP

SUBROUTINE GetAFromM(M, AStar, A)
  REAL(KIND = rDef), INTENT(IN):: M, AStar
  REAL(KIND = rDef), INTENT(OUT):: A

! local variables

  REAL(KIND = rDef):: fac, pow, gm1, gp1, gm1o2, aOverAStar2

  gp1 = gam+1.0_rDef
  gm1 = gam-1.0_rDef
  gm1o2 = 0.5_rDef*gm1

  pow = gp1/gm1

  fac = (2.0_rDef/gp1)*(1.0_rDef+gm1o2*m*m)
  aOverAStar2 = (1.0_rDef/(M*M)) &
                *(fac**pow) 
  
  A = AStar*SQRT(aOverAStar2)

  RETURN 
END SUBROUTINE GetAFromM

SUBROUTINE NormalShockRelationsv1(M1,  &
                                  p01, &
                                  t01, &
                                  M2,  &
                                  p02, &
                                  t02)

  REAL(KIND = rDef), INTENT(IN)  :: M1, p01, t01
  REAL(KIND = rDef), INTENT(OUT):: M2, p02, t02

! local data

  REAL(KIND = rDef):: num1, den1, gm1, pow1, &
                     num2, den2, gp1, pow2

  gm1 = gam-1.0_rDef
  gp1 = gam+1.0_rDef

  t02 = t01

  num1 = (gm1*M1*M1) + 2.0_rDef
  den1 = (2.0_rDef*gam*M1*M1) - gm1
  M2   = SQRT(num1/den1)

  pow1 = gam/gm1
  pow2 = 1.0_rDef/gm1

  den2 = den1
  den1 = num1
  num1 = gp1*M1*M1
  num2 = gp1

  p02 = p01*((num1/den1)**pow1)*((num2/den2)**pow2)

  RETURN
END SUBROUTINE NormalShockRelationsv1

END MODULE Compressible1DFunctions
