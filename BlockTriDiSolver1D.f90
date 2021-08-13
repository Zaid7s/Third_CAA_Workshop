MODULE BlockTriDiSolver1D

    USE,  INTRINSIC :: ISO_FORTRAN_ENV
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: BlockTriMatrix

    INTEGER, PARAMETER :: rDef = REAL64

    CONTAINS
 
    SUBROUTINE BlockTriMatrix(  numVar,         &
                                iSolveStart,    &
                                iSolveEnd,      &
                                aStart,         &
                                A,              &
                                bStart,         &
                                B,              &
                                cStart,         &
                                C,              &
                                xStart,         &
                                X,              &
                                rStart,         &
                                R)


    INTEGER, INTENT(IN) :: numVar
    INTEGER, DIMENSION(:), INTENT(IN) :: iSolveStart,   & ! dimension 1
                                           iSolveEnd,   &
                                           aStart,      &
                                           bStart,      &
                                           cStart,      &
                                           xStart,      &
                                           rStart
    REAL(KIND=rDef), DIMENSION(aStart(1):,:,:), INTENT(INOUT) :: A
    REAL(KIND=rDef), DIMENSION(bStart(1):,:,:), INTENT(INOUT) :: B
    REAL(KIND=rDef), DIMENSION(cStart(1):,:,:), INTENT(INOUT) :: C
    REAL(KIND=rDef), DIMENSION(xStart(1):,:), INTENT(INOUT) :: X
    REAL(KIND=rDef), DIMENSION(rStart(1):,:), INTENT(INOUT) :: R

    INTEGER :: I, M, N, NN
    REAL(KIND=rDef) :: S, T


    CONTINUE ! execution begins here

        DO I=iSolveStart(1),iSolveEnd(1)-1
            DO M=1,numVar-1
                T=1.0_rDef/B(I,M,M)
                
                DO N=M+1,numVar
                    B(I,M,N) = T*B(I,M,N)
                END DO

                DO N=1,numVar
                    C(I,M,N) = T*C(I,M,N)
                END DO

                R(I,M) = T*R(I,M)

                DO N=M+1,numVar
                    DO NN=M+1,numVar
                        B(I,N,NN) = B(I,N,NN)-B(I,N,M)*B(I,M,NN)
                    END DO

                    DO NN=1,numVar
                        C(I,N,NN) = C(I,N,NN)-B(I,N,M)*C(I,M,NN)
                    END DO
                    R(I,N) = R(I,N)-B(I,N,M)*R(I,M)
                END DO ! n loop
            END DO ! m loop

            T= 1.0_rDef/B(I,numVar,numVar)

            DO NN = 1,numVar
                C(I,numVar,NN) = T*C(I,numVar,NN)
            END DO

            R(I,numVar) = T*R(I,numVar)

            DO N=numVar,2,-1
                DO M=N-1,1,-1
                    DO NN=1,numVar
                        C(I,M,NN) = C(I,M,NN)-B(I,M,N)*C(I,N,NN)
                    END DO

                    R(I,M) = R(I,M)-B(I,M,N)*R(I,N)
                END DO ! m loop
            END DO  ! n loop

            DO N=1,numVar
                DO M=1,numVar
                    DO NN=1,numVar
                        B(I+1,M,NN) = B(I+1,M,NN)-A(I+1,M,N)*C(I,N,NN)
                    END DO ! nn loop
                    R(I+1,M) = R(I+1,M)-A(I+1,M,N)*R(I,N)
                END DO ! m loop
            END DO ! n loop
        END DO  ! I loop


        I = iSolveEnd(1)


        DO M=1,numVar-1

            T=1.0_rDef/B(I,M,M)

            DO N=M+1,numVar
                B(I,M,N) = T*B(I,M,N)
            END DO

            R(I,M) = T*R(I,M)


            DO N=M+1,numVar

                DO NN=M+1,numVar
                    B(I,N,NN) = B(I,N,NN)-B(I,N,M)*B(I,M,NN)
                END DO

                R(I,N) = R(I,N)-B(I,N,M)*R(I,M)

            END DO ! n loop

        END DO ! m loop

        T= 1.0_rDef/B(I,numVar,numVar)

        R(I,numVar) = T*R(I,numVar)

        DO N=numVar,2,-1
            DO M=N-1,1,-1
                R(I,M) = R(I,M)-B(I,M,N)*R(I,N)
            END DO ! m loop
        END DO  ! n loop

        DO M=1,numVar
            X(I,M) = R(I,M)
        END DO

        DO I=iSolveEnd(1)-1,iSolveStart(1),-1
            DO M=numVar,1,-1
                S=0.0_rDef
                DO NN=1,numVar
                    S = C(I,M,NN)*X(I+1,NN)+S
                END DO ! nn loop
                X(I,M) = R(I,M)-S
            END DO ! m loop
        END DO ! i loop

      RETURN

    END SUBROUTINE BlockTriMatrix

END MODULE BlockTriDiSolver1D
